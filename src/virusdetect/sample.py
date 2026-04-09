from __future__ import annotations

import gzip
from pathlib import Path


def open_maybe_gzip(path: Path, mode: str = "rt"):
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8")
    return path.open(mode, encoding="utf-8")


def detect_file_type(path: str | Path) -> str:
    candidate = Path(path)
    with open_maybe_gzip(candidate) as handle:
        first_line = handle.readline()

    if first_line.startswith(">"):
        return "fasta"
    if first_line.startswith("@"):
        return "fastq"
    raise ValueError(f"Cannot detect file type for {candidate}")


def iter_sequences(path: str | Path):
    candidate = Path(path)
    file_type = detect_file_type(candidate)
    with open_maybe_gzip(candidate) as handle:
        while True:
            header = handle.readline()
            if not header:
                break

            sequence = handle.readline()
            if not sequence:
                raise ValueError(f"Truncated sequence record in {candidate}")

            header = header.rstrip("\n")
            sequence = sequence.rstrip("\n")

            if file_type == "fastq":
                plus = handle.readline()
                quality = handle.readline()
                if not plus or not quality:
                    raise ValueError(f"Truncated FASTQ record in {candidate}")
                yield header, sequence, plus, quality
            else:
                yield header, sequence, None, None


def detect_data_type(path: str | Path) -> str:
    total_length = 0
    count = 0
    for _, sequence, _, _ in iter_sequences(path):
        total_length += len(sequence)
        count += 1
        if count >= 1000:
            break

    if count == 0:
        raise ValueError(f"No sequences found in {path}")

    mean_length = total_length / count
    return "mRNA" if mean_length >= 36 else "sRNA"


def count_sequences(path: str | Path) -> int:
    count = 0
    for _record in iter_sequences(path):
        count += 1
    return count


def parse_read_length_spec(spec: str) -> set[int]:
    selected: set[int] = set()
    for segment in spec.split(":"):
        segment = segment.strip()
        if not segment:
            continue
        if "-" in segment:
            start_text, end_text = segment.split("-", 1)
            start = int(start_text)
            end = int(end_text)
            if start > end:
                raise ValueError(f"Invalid read-length range: {segment}")
            for length in range(start, end + 1):
                selected.add(length)
        else:
            selected.add(int(segment))

    if not selected or min(selected) <= 0:
        raise ValueError(f"Invalid read-length spec: {spec}")
    return selected


def filter_reads_by_length(read_length: str, input_path: str | Path, output_path: str | Path) -> int:
    selected = parse_read_length_spec(read_length)
    input_candidate = Path(input_path)
    output_candidate = Path(output_path)
    file_type = detect_file_type(input_candidate)
    kept = 0

    with output_candidate.open("w", encoding="utf-8") as output_handle:
        for header, sequence, plus, quality in iter_sequences(input_candidate):
            if len(sequence) not in selected:
                continue
            output_handle.write(f"{header}\n{sequence}\n")
            if file_type == "fastq":
                output_handle.write(plus)
                output_handle.write(quality)
            kept += 1

    return kept

