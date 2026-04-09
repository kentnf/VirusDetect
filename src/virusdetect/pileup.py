from __future__ import annotations

import gzip
from pathlib import Path


def open_text(path: str | Path):
    candidate = Path(path)
    if candidate.suffix == ".gz":
        return gzip.open(candidate, "rt", encoding="utf-8")
    return candidate.open("r", encoding="utf-8")


def load_seq_lengths_from_info(seq_info_path: str | Path, ref_ids: set[str]) -> dict[str, int]:
    seq_lengths: dict[str, int] = {ref_id: 0 for ref_id in ref_ids}
    with open_text(seq_info_path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if fields[0] in seq_lengths:
                seq_lengths[fields[0]] = int(fields[1])
    return seq_lengths


def filter_pileup_by_coverage(input_pileup: str | Path, seq_info_path: str | Path, coverage: float, output_pileup: str | Path) -> None:
    ref_ids: set[str] = set()
    with Path(input_pileup).open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) >= 6:
                ref_ids.add(fields[0])

    seq_lengths = load_seq_lengths_from_info(seq_info_path, ref_ids)

    pre_id = None
    pileup_info: list[str] = []
    with Path(input_pileup).open("r", encoding="utf-8") as input_handle, Path(output_pileup).open("w", encoding="utf-8") as output_handle:
        for raw_line in input_handle:
            line = raw_line.rstrip("\n")
            fields = line.split("\t")
            if len(fields) < 6:
                continue
            ref_id = fields[0]

            if pre_id is not None and ref_id != pre_id:
                len_cutoff = seq_lengths[pre_id] * coverage
                if len(pileup_info) > len_cutoff:
                    output_handle.write("\n".join(pileup_info) + "\n")
                pileup_info = []

            pileup_info.append(line)
            pre_id = ref_id

        if pre_id is not None:
            len_cutoff = seq_lengths[pre_id] * coverage
            if len(pileup_info) > len_cutoff:
                output_handle.write("\n".join(pileup_info) + "\n")


def load_fasta_lengths(fasta_path: str | Path) -> dict[str, int]:
    seq_lengths: dict[str, int] = {}
    current_id = None
    current_length = 0
    with Path(fasta_path).open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    seq_lengths[current_id] = current_length
                current_id = line[1:].split()[0]
                current_length = 0
            else:
                current_length += len(line.strip())
        if current_id is not None:
            seq_lengths[current_id] = current_length
    return seq_lengths


def parse_indel_str(map_str: str, i: int) -> tuple[int, str]:
    i += 1
    sub_map_str = map_str[i:]
    ins_num_chars: list[str] = []
    for char in sub_map_str:
        if char.isdigit():
            ins_num_chars.append(char)
        else:
            break
    ins_num = int("".join(ins_num_chars))
    begin = i + len(ins_num_chars)
    end = i + ins_num
    ins_str = map_str[begin:begin + ins_num]
    return end, ins_str


def get_best_base(ref_base: str, depth: int, align: str) -> str:
    acgt = ("A", "C", "G", "T")
    match_num = 0
    ins_num = 0
    del_num = 0
    depth_float = float(depth) + 0.00001
    mut_depth = depth_float
    ins_depth = depth_float
    del_depth = depth_float

    ref_base = ref_base.upper()
    best_base = ref_base
    map_str = align.upper()
    base_count = {base: 0 for base in acgt}
    ins_list: list[str] = []
    i = 0
    while i < len(map_str):
        letter = map_str[i]
        if letter == "^":
            i += 1
            del_depth -= 1
        elif letter == "$":
            ins_depth -= 1
            del_depth -= 1
        elif letter == "+":
            i, ins_str = parse_indel_str(map_str, i)
            ins_list.append(ins_str)
            ins_num += 1
        elif letter == "-":
            i, _ins_str = parse_indel_str(map_str, i)
        elif letter == "*":
            del_num += 1
        elif letter in {",", "."}:
            match_num += 1
        elif letter in {"N", "n"}:
            mut_depth -= 1
        elif letter in base_count:
            base_count[letter] += 1
        i += 1

    base_count[ref_base] = match_num
    max_mut_count = 0
    for base in acgt:
        if base_count[base] > max_mut_count:
            max_mut_count = base_count[base]
            best_base = base

    max_mut_percent = max_mut_count / mut_depth if mut_depth > 0 else 0
    ins_percent = ins_num / ins_depth if ins_depth > 0 else 0
    del_percent = del_num / del_depth if del_depth > 0 else 0

    max_count = 0
    max_str = ""
    ins_count: dict[str, int] = {}
    for ins_str in ins_list:
        ins_count[ins_str] = ins_count.get(ins_str, 0) + 1
    for ins_str in ins_list:
        if ins_count[ins_str] > max_count:
            max_count = ins_count[ins_str]
            max_str = ins_str

    max_percent = max_mut_percent
    if del_percent > max_percent:
        max_percent = del_percent
        best_base = ""
    if ins_percent >= max_percent:
        best_base = ref_base + max_str

    return best_base


def pileup_to_contig(
    input_seq: str | Path,
    input_pileup: str | Path,
    output_seq: str | Path,
    len_cutoff: int,
    depth_cutoff: float,
    prefix: str = "PILEUP2SEQ",
) -> int:
    seq_lengths = load_fasta_lengths(input_seq)
    input_path = Path(input_pileup)
    output_path = Path(output_seq)

    first_record = None
    with input_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            fields = line.split("\t")
            if len(fields) < 6:
                continue
            ref_id, pos, ref_base, depth, align = fields[0], int(fields[1]), fields[2], int(fields[3]), fields[4]
            first_record = (ref_id, pos, get_best_base(ref_base, depth, align), depth)
            break

    if first_record is None:
        output_path.write_text("", encoding="utf-8")
        return 0

    contig_num = 0
    pre_id, pre_pos, first_base, dep = first_record
    seq = first_base
    length = 1

    with input_path.open("r", encoding="utf-8") as input_handle, output_path.open("w", encoding="utf-8") as output_handle:
        first_consumed = False
        for raw_line in input_handle:
            line = raw_line.rstrip("\n")
            fields = line.split("\t")
            if len(fields) < 6:
                continue
            ref_id, pos, ref_base, depth, align = fields[0], int(fields[1]), fields[2], int(fields[3]), fields[4]

            if not first_consumed:
                first_consumed = True
                continue

            if pos > seq_lengths.get(ref_id, 0):
                continue

            best_base = get_best_base(ref_base, depth, align)
            if ref_id == pre_id:
                if pos == pre_pos + 1:
                    length += 1
                    dep += depth
                    seq += best_base
                    pre_pos = pos
                else:
                    if length > len_cutoff and (dep / length >= depth_cutoff):
                        contig_num += 1
                        output_handle.write(f">{prefix}{contig_num} {pre_id} {pre_pos}\n{seq}\n")
                    length = 1
                    dep = depth
                    seq = best_base
                    pre_pos = pos
                    pre_id = ref_id
            else:
                if length > len_cutoff and (dep / length >= depth_cutoff):
                    contig_num += 1
                    output_handle.write(f">{prefix}{contig_num} {pre_id} {pre_pos}\n{seq}\n")
                length = 1
                dep = depth
                seq = best_base
                pre_pos = pos
                pre_id = ref_id

        if length > len_cutoff and (dep / length > depth_cutoff):
            contig_num += 1
            output_handle.write(f">{prefix}{contig_num} {pre_id} {pre_pos}\n{seq}\n")

    return contig_num


def count_fasta_sequences(fasta_path: str | Path) -> int:
    count = 0
    with Path(fasta_path).open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            if raw_line.startswith(">"):
                count += 1
    return count
