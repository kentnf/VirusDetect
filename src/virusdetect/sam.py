from __future__ import annotations

import re
from pathlib import Path

NM_RE = re.compile(r"\tNM:i:(\d+)")
XA_RE = re.compile(r"\tXA:Z:(\S+)")


def reverse_complement(sequence: str) -> str:
    table = str.maketrans("ACGTacgt", "TGCAtgca")
    return sequence.translate(table)[::-1]


def expand_xa_hits(input_sam: str | Path) -> None:
    input_path = Path(input_sam)
    temp_path = Path(f"{input_path}.temp")

    with input_path.open("r", encoding="utf-8") as input_handle, temp_path.open("w", encoding="utf-8") as output_handle:
        for raw_line in input_handle:
            xa_match = XA_RE.search(raw_line)
            if xa_match is None:
                output_handle.write(raw_line)
                continue

            output_handle.write(raw_line)
            fields = raw_line.rstrip("\n").split("\t")
            xa_blob = xa_match.group(1)

            for alt_ref, alt_pos, cigar, nm_text in re.findall(r"([^,;]+),([-+]\d+),([^,]+),(\d+);", xa_blob):
                sequence = fields[9]
                quality = fields[10]
                is_primary_reverse = (int(fields[1]) & 0x10) > 0
                is_alt_reverse = int(alt_pos) < 0

                if is_primary_reverse ^ is_alt_reverse:
                    sequence = reverse_complement(sequence)
                    quality = quality[::-1]

                alt_flag = (int(fields[1]) & 0x6E9) | (0x10 if is_alt_reverse else 0)
                alt_record = [
                    fields[0],
                    str(alt_flag),
                    alt_ref,
                    str(abs(int(alt_pos))),
                    "0",
                    cigar,
                    *fields[6:8],
                    "0",
                    sequence,
                    quality,
                    f"NM:i:{nm_text}",
                ]
                output_handle.write("\t".join(alt_record) + "\n")

    temp_path.replace(input_path)


def edit_distance_from_sam(line: str) -> int:
    match = NM_RE.search(line)
    if match is None:
        raise ValueError(f"SAM alignment is missing NM tag: {line.rstrip()}")
    return int(match.group(1))


def filter_sam_best_hits(input_sam: str | Path, max_distance: int = 2) -> int:
    input_path = Path(input_sam)
    temp_path = Path(f"{input_path}.temp")

    mapped_records: list[str] = []
    with input_path.open("r", encoding="utf-8") as input_handle:
        for raw_line in input_handle:
            if raw_line.startswith("@"):
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            if fields[1] == "4":
                continue
            mapped_records.append(raw_line.rstrip("\n"))

    kept_records: list[str] = []
    current_query = None
    current_alignments: list[str] = []
    best_distance = -1

    def flush_current():
        nonlocal current_alignments, best_distance
        if not current_alignments:
            return
        for alignment in current_alignments:
            if edit_distance_from_sam(alignment) == best_distance:
                kept_records.append(alignment)
        current_alignments = []
        best_distance = -1

    for alignment in mapped_records:
        fields = alignment.split("\t")
        query_name = fields[0]
        if query_name != current_query:
            flush_current()
            current_query = query_name

        distance = edit_distance_from_sam(alignment)
        if distance >= max_distance:
            continue
        if best_distance == -1 or distance < best_distance:
            best_distance = distance
        current_alignments.append(alignment)

    flush_current()

    with temp_path.open("w", encoding="utf-8") as output_handle:
        for alignment in kept_records:
            output_handle.write(alignment + "\n")

    temp_path.replace(input_path)
    return len({line.split("\t", 1)[0] for line in kept_records})


def generate_unmapped_reads_from_sam(input_sam: str | Path, output_reads: str | Path) -> int:
    input_path = Path(input_sam)
    output_path = Path(output_reads)
    unmapped_num = 0

    with input_path.open("r", encoding="utf-8") as input_handle, output_path.open("w", encoding="utf-8") as output_handle:
        for raw_line in input_handle:
            if raw_line.startswith("@"):
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            if fields[1] == "4":
                output_handle.write(f">{fields[0]}\n{fields[9]}\n")
                unmapped_num += 1

    return unmapped_num
