from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class FastaRecord:
    header: str
    sequence: str

    @property
    def seq_id(self) -> str:
        return self.header.split()[0]


def read_fasta_records(path: str | Path) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    header: str | None = None
    sequence_parts: list[str] = []

    with Path(path).open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append(FastaRecord(header=header, sequence="".join(sequence_parts)))
                header = line[1:].strip()
                sequence_parts = []
            else:
                sequence_parts.append(line.strip())

    if header is not None:
        records.append(FastaRecord(header=header, sequence="".join(sequence_parts)))

    return records


def write_fasta_records(path: str | Path, records: list[FastaRecord]) -> None:
    output_path = Path(path)
    with output_path.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(f">{record.header}\n{record.sequence}\n")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTacgtNn", "TGCAtgcaNn"))[::-1]
