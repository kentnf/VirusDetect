from __future__ import annotations

import csv
import gzip
import html
import json
import shutil
import subprocess
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path

from virusdetect.db import resolve_data_path, resolve_reference_path, resolve_seq_info_path
from virusdetect.fasta import FastaRecord, read_fasta_records, write_fasta_records
from virusdetect.runtime import tool_path_map
from virusdetect.sample import count_sequences
from virusdetect.sam import expand_xa_hits, reverse_complement

BLAST_OUTFMT_FIELDS = (
    "qseqid",
    "qlen",
    "sseqid",
    "slen",
    "length",
    "pident",
    "evalue",
    "bitscore",
    "qstart",
    "qend",
    "sstart",
    "send",
    "qcovs",
    "qseq",
    "sseq",
)
BLASTN_MIN_IDENTITY = 60.0
BLASTN_MIN_QUERY_COVERAGE = 50.0
BLASTX_MIN_QUERY_COVERAGE = 60.0
ALIGNMENT_WRAP_WIDTH = 95


@dataclass(frozen=True)
class VirusAnnotation:
    reference_id: str
    hit_length: int
    genus: str
    description: str


@dataclass(frozen=True)
class BlastHit:
    analysis: str
    contig_id: str
    contig_length: int
    hit_id: str
    reference_id: str
    hit_length: int
    genus: str
    description: str
    alignment_length: int
    percent_identity: float
    evalue: float
    bit_score: float
    query_start: int
    query_end: int
    hit_start: int
    hit_end: int
    query_coverage: float
    query_sequence: str = ""
    hit_sequence: str = ""


@dataclass(frozen=True)
class ContigDepth:
    contig_id: str
    contig_length: int
    covered_bases: int
    total_depth: int
    mean_depth: float
    normalized_depth: float


@dataclass(frozen=True)
class IdentifyResult:
    result_dir: str
    blastn_raw_tsv: str
    blastn_top_tsv: str
    blastx_raw_tsv: str
    blastx_top_tsv: str
    summary_tsv: str
    summary_json: str
    known_contig_count: int
    novel_contig_count: int
    undetermined_contig_count: int


def run_command(command: list[str], debug: bool = False) -> None:
    if debug:
        print("CMD:", " ".join(command))
    completed = subprocess.run(command, check=False)
    if completed.returncode != 0:
        raise SystemExit(f"Command failed with exit code {completed.returncode}: {' '.join(command)}")


def run_command_to_file(command: list[str], output_path: Path, binary: bool = False, debug: bool = False) -> None:
    if debug:
        print("CMD:", " ".join(command), ">", output_path)
    mode = "wb" if binary else "w"
    kwargs = {} if binary else {"encoding": "utf-8"}
    with output_path.open(mode, **kwargs) as handle:
        completed = subprocess.run(command, stdout=handle, check=False)
    if completed.returncode != 0:
        raise SystemExit(f"Command failed with exit code {completed.returncode}: {' '.join(command)}")


def open_text_maybe_gzip(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def normalize_seq_id(raw_value: str) -> str:
    token = raw_value.strip().split()[0]
    if token.startswith("lcl|"):
        token = token[4:]
    if "|" not in token:
        return token

    parts = [part for part in token.split("|") if part]
    for marker in ("ref", "gb", "emb", "dbj", "sp", "tr"):
        if marker in parts:
            index = parts.index(marker)
            if index + 1 < len(parts):
                return parts[index + 1]
    return parts[-1]


def load_database_annotations(db_path: str) -> tuple[dict[str, VirusAnnotation], dict[str, str]]:
    annotations: dict[str, VirusAnnotation] = {}
    seq_info_path = Path(resolve_seq_info_path(db_path))
    with open_text_maybe_gzip(seq_info_path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                continue
            reference_id = normalize_seq_id(fields[0])
            try:
                hit_length = int(fields[1])
            except ValueError:
                hit_length = 0
            annotations[reference_id] = VirusAnnotation(
                reference_id=reference_id,
                hit_length=hit_length,
                genus=fields[2],
                description=fields[3],
            )

    protein_to_reference: dict[str, str] = {}
    mapping_path = Path(resolve_data_path("vrl_idmapping.gz", db_path))
    with open_text_maybe_gzip(mapping_path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                continue
            protein_to_reference[normalize_seq_id(fields[1])] = normalize_seq_id(fields[0])

    return annotations, protein_to_reference


def has_blast_db(prefix: Path, dbtype: str) -> bool:
    if dbtype == "nucl":
        suffixes = (".nhr", ".nin", ".nsq", ".ndb", ".00.nin")
    else:
        suffixes = (".phr", ".pin", ".psq", ".pdb", ".00.pin")
    return any((prefix.parent / f"{prefix.name}{suffix}").exists() for suffix in suffixes)


def ensure_blast_db(reference_path: str, dbtype: str, tool_paths: dict[str, str], debug: bool = False) -> None:
    reference_prefix = Path(reference_path)
    if has_blast_db(reference_prefix, dbtype):
        return

    makeblastdb = tool_paths.get("makeblastdb")
    if not makeblastdb:
        raise SystemExit("makeblastdb is required for the Python identify stage")

    run_command(
        [
            makeblastdb,
            "-in",
            str(reference_prefix),
            "-out",
            str(reference_prefix),
            "-dbtype",
            dbtype,
            "-parse_seqids",
        ],
        debug=debug,
    )


def ensure_bwa_index(reference_path: str, tool_paths: dict[str, str], debug: bool = False) -> None:
    if Path(f"{reference_path}.amb").exists():
        return

    bwa = tool_paths.get("bwa")
    if not bwa:
        raise SystemExit("bwa is required for contig-depth estimation")

    run_command([bwa, "index", "-p", reference_path, "-a", "bwtsw", reference_path], debug=debug)


def ensure_fasta_index(reference_path: str, tool_paths: dict[str, str], debug: bool = False) -> None:
    if Path(f"{reference_path}.fai").exists():
        return

    samtools = tool_paths.get("samtools")
    if not samtools:
        raise SystemExit("samtools is required for contig-depth estimation")

    run_command([samtools, "faidx", reference_path], debug=debug)


def resolve_protein_reference_path(reference: str, db_path: str) -> str:
    explicit_candidate = Path(f"{reference}_prot")
    if explicit_candidate.exists():
        return str(explicit_candidate.resolve())

    reference_name = Path(reference).name
    return resolve_data_path(f"{reference_name}_prot", db_path)


def build_bwa_align_command(reference_path: str, sample_path: str, threads: int) -> list[str]:
    return [
        "bwa",
        "aln",
        "-n",
        "1",
        "-o",
        "1",
        "-e",
        "1",
        "-i",
        "0",
        "-l",
        "15",
        "-k",
        "1",
        "-t",
        str(threads),
        reference_path,
        sample_path,
    ]


def build_bwa_samse_command(reference_path: str, sai_path: str, sample_path: str) -> list[str]:
    return [
        "bwa",
        "samse",
        "-n",
        "10000",
        reference_path,
        sai_path,
        sample_path,
    ]


def collect_contig_read_length_stats(sam_path: Path) -> dict[str, dict[int, int]]:
    stats: dict[str, dict[int, int]] = defaultdict(dict)
    with sam_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            if not raw_line or raw_line.startswith("@"):
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            try:
                flag = int(fields[1])
            except ValueError:
                continue
            contig_id = fields[2]
            sequence = fields[9]
            if contig_id == "*" or flag & 0x4 or sequence == "*":
                continue
            read_length = len(sequence)
            stats.setdefault(contig_id, {})
            stats[contig_id][read_length] = stats[contig_id].get(read_length, 0) + 1
    return {contig_id: dict(sorted(length_counts.items())) for contig_id, length_counts in stats.items()}


def cleanup_intermediate_files(paths: list[Path], keep_files: bool = False) -> None:
    if keep_files:
        return

    for path in paths:
        try:
            path.unlink()
        except FileNotFoundError:
            continue


def cleanup_result_dir(output_dir: Path, keep_names: set[str]) -> None:
    for path in output_dir.iterdir():
        if path.name in keep_names:
            continue
        if path.is_dir():
            shutil.rmtree(path, ignore_errors=True)
        else:
            try:
                path.unlink()
            except FileNotFoundError:
                continue


def result_file(path_dir: Path, name: str) -> Path:
    return path_dir / name


def compute_contig_depths(
    sample_path: Path,
    contig_path: Path,
    contig_records: list[FastaRecord],
    result_dir: Path,
    args,
    tool_paths: dict[str, str],
) -> tuple[dict[str, ContigDepth], Path, dict[str, dict[int, int]]]:
    depth_path = result_file(result_dir, "contig_depth.tsv")
    if not contig_records:
        write_contig_depth_tsv(depth_path, {})
        return {}, depth_path, {}

    bwa = tool_paths.get("bwa")
    samtools = tool_paths.get("samtools")
    if not bwa or not samtools:
        raise SystemExit("bwa and samtools are required for contig-depth estimation")

    reference_path = str(contig_path)
    ensure_bwa_index(reference_path, {"bwa": bwa}, debug=args.debug)
    ensure_fasta_index(reference_path, {"samtools": samtools}, debug=args.debug)

    sai_path = result_file(result_dir, "contig_depth.bwa.sai")
    sam_path = result_file(result_dir, "contig_depth.sam")
    bam_path = result_file(result_dir, "contig_depth.bam")
    sorted_prefix = result_file(result_dir, "contig_depth.sorted")
    sorted_bam_path = result_file(result_dir, "contig_depth.sorted.bam")
    raw_depth_path = result_file(result_dir, "contig_depth.pileup")

    align_command = build_bwa_align_command(reference_path, str(sample_path), args.threads)
    align_command[0] = bwa
    run_command_to_file(align_command, sai_path, debug=args.debug)

    samse_command = build_bwa_samse_command(reference_path, str(sai_path), str(sample_path))
    samse_command[0] = bwa
    run_command_to_file(samse_command, sam_path, debug=args.debug)
    expand_xa_hits(sam_path)
    read_length_stats = collect_contig_read_length_stats(sam_path)

    view_command = [samtools, "view", "-bt", f"{reference_path}.fai", str(sam_path)]
    run_command_to_file(view_command, bam_path, binary=True, debug=args.debug)
    run_command([samtools, "sort", str(bam_path), str(sorted_prefix)], debug=args.debug)
    if not sorted_bam_path.exists():
        raise SystemExit(f"samtools sort did not produce expected output: {sorted_bam_path}")
    run_command_to_file([samtools, "mpileup", "-f", reference_path, str(sorted_bam_path)], raw_depth_path, debug=args.debug)

    contig_lengths = {record.seq_id: len(record.sequence) for record in contig_records}
    covered_bases = {contig_id: 0 for contig_id in contig_lengths}
    total_depth = {contig_id: 0 for contig_id in contig_lengths}

    with raw_depth_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 4 or fields[0] not in contig_lengths:
                continue
            depth = int(fields[3])
            total_depth[fields[0]] += depth
            if depth > 0:
                covered_bases[fields[0]] += 1

    library_size = count_sequences(sample_path)
    contig_depths: dict[str, ContigDepth] = {}
    for contig_id, contig_length in sorted(contig_lengths.items()):
        total = total_depth[contig_id]
        mean_depth = (total / contig_length) if contig_length else 0.0
        normalized_depth = ((1e6 * total) / (contig_length * library_size)) if contig_length and library_size else 0.0
        contig_depths[contig_id] = ContigDepth(
            contig_id=contig_id,
            contig_length=contig_length,
            covered_bases=covered_bases[contig_id],
            total_depth=total,
            mean_depth=mean_depth,
            normalized_depth=normalized_depth,
        )

    write_contig_depth_tsv(depth_path, contig_depths)
    cleanup_intermediate_files(
        [
            sai_path,
            sam_path,
            bam_path,
            sorted_bam_path,
            raw_depth_path,
        ],
        keep_files=args.debug,
    )
    return contig_depths, depth_path, read_length_stats


def run_blastn(
    query_path: Path,
    reference_path: str,
    output_path: Path,
    args,
    tool_paths: dict[str, str],
) -> None:
    blastn = tool_paths.get("blastn")
    if not blastn:
        raise SystemExit("blastn is required for the Python identify stage")

    command = [
        blastn,
        "-task",
        "megablast",
        "-query",
        str(query_path),
        "-db",
        reference_path,
        "-out",
        str(output_path),
        "-outfmt",
        "6 " + " ".join(BLAST_OUTFMT_FIELDS),
        "-word_size",
        str(args.word_size),
        "-evalue",
        str(args.exp_value),
        "-dust",
        "no",
        "-max_target_seqs",
        "25",
        "-num_threads",
        str(args.threads),
    ]
    run_command(command, debug=args.debug)


def run_blastx(
    query_path: Path,
    reference_path: str,
    output_path: Path,
    args,
    effective_exp_valuex: float,
    tool_paths: dict[str, str],
) -> None:
    blastx = tool_paths.get("blastx")
    if not blastx:
        raise SystemExit("blastx is required for the Python identify stage")

    command = [
        blastx,
        "-query",
        str(query_path),
        "-db",
        reference_path,
        "-out",
        str(output_path),
        "-outfmt",
        "6 " + " ".join(BLAST_OUTFMT_FIELDS),
        "-evalue",
        str(effective_exp_valuex),
        "-seg",
        "no",
        "-max_target_seqs",
        "25",
        "-num_threads",
        str(args.threads),
    ]
    run_command(command, debug=args.debug)


def parse_blast_hits(
    analysis: str,
    output_path: Path,
    annotations: dict[str, VirusAnnotation],
    protein_to_reference: dict[str, str],
) -> list[BlastHit]:
    hits: list[BlastHit] = []
    if not output_path.exists() or output_path.stat().st_size == 0:
        return hits

    with output_path.open("r", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            hit_id = normalize_seq_id(row[2])
            reference_id = hit_id
            if analysis == "blastx":
                reference_id = protein_to_reference.get(hit_id, hit_id)
            annotation = annotations.get(reference_id)
            hit_length = int(row[3])
            hits.append(
                BlastHit(
                    analysis=analysis,
                    contig_id=row[0],
                    contig_length=int(row[1]),
                    hit_id=hit_id,
                    reference_id=reference_id,
                    hit_length=hit_length,
                    genus=annotation.genus if annotation else "",
                    description=annotation.description if annotation else "",
                    alignment_length=int(row[4]),
                    percent_identity=float(row[5]),
                    evalue=float(row[6]),
                    bit_score=float(row[7]),
                    query_start=int(row[8]),
                    query_end=int(row[9]),
                    hit_start=int(row[10]),
                    hit_end=int(row[11]),
                    query_coverage=float(row[12]),
                    query_sequence=row[13] if len(row) > 13 else "",
                    hit_sequence=row[14] if len(row) > 14 else "",
                )
            )

    return hits


def hit_sort_key(hit: BlastHit) -> tuple[float, float, float, float, int, str]:
    return (
        hit.evalue,
        -hit.bit_score,
        -hit.query_coverage,
        -hit.percent_identity,
        -hit.alignment_length,
        hit.hit_id,
    )


def hit_group_id(hit: BlastHit) -> str:
    return hit.reference_id if hit.analysis == "blastn" else hit.hit_id


def merge_covered_bases(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0

    merged_total = 0
    current_start, current_end = sorted((min(intervals[0]), max(intervals[0])))
    for start, end in intervals[1:]:
        start, end = sorted((start, end))
        if start <= current_end + 1:
            if end > current_end:
                current_end = end
            continue
        merged_total += current_end - current_start + 1
        current_start, current_end = start, end

    merged_total += current_end - current_start + 1
    return merged_total


def filter_hits(hits: list[BlastHit], min_identity: float, min_query_coverage: float) -> list[BlastHit]:
    return [
        hit
        for hit in hits
        if hit.percent_identity >= min_identity and hit.query_coverage >= min_query_coverage
    ]


def group_hits_by_reference(hits: list[BlastHit]) -> dict[tuple[str, str, str, int, str, str], list[BlastHit]]:
    grouped: dict[tuple[str, str, str, int, str, str], list[BlastHit]] = defaultdict(list)
    for hit in hits:
        grouped[(hit.analysis, hit_group_id(hit), hit.reference_id, hit.hit_length, hit.genus, hit.description)].append(hit)
    return grouped


def select_best_hits_by_contig(hits: list[BlastHit]) -> dict[str, BlastHit]:
    best_by_contig: dict[str, BlastHit] = {}
    for hit in hits:
        current = best_by_contig.get(hit.contig_id)
        if current is None or hit_sort_key(hit) < hit_sort_key(current):
            best_by_contig[hit.contig_id] = hit
    return best_by_contig


def select_best_hits(
    hits: list[BlastHit],
    min_identity: float,
    min_query_coverage: float,
    allowed_group_ids: set[str] | None = None,
) -> list[BlastHit]:
    selected: dict[str, BlastHit] = {}
    for hit in filter_hits(hits, min_identity, min_query_coverage):
        if allowed_group_ids is not None and hit_group_id(hit) not in allowed_group_ids:
            continue
        current = selected.get(hit.contig_id)
        if current is None or hit_sort_key(hit) < hit_sort_key(current):
            selected[hit.contig_id] = hit
    return [selected[contig_id] for contig_id in sorted(selected)]


def write_hit_table(path: Path, hits: list[BlastHit]) -> None:
    field_names = [
        "analysis",
        "contig_id",
        "contig_length",
        "hit_id",
        "reference_id",
        "hit_length",
        "genus",
        "description",
        "alignment_length",
        "percent_identity",
        "evalue",
        "bit_score",
        "query_start",
        "query_end",
        "hit_start",
        "hit_end",
        "query_coverage",
        "query_sequence",
        "hit_sequence",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, field_names, delimiter="\t")
        writer.writeheader()
        for hit in hits:
            writer.writerow(asdict(hit))


def compute_alignment_identity_text(hit: BlastHit) -> str:
    query_sequence = hit.query_sequence.replace(" ", "")
    hit_sequence = hit.hit_sequence.replace(" ", "")
    if not query_sequence or not hit_sequence:
        matched = int(round(hit.alignment_length * (hit.percent_identity / 100.0)))
        return f"{matched}/{hit.alignment_length}({hit.percent_identity:.0f}%)"

    matched = 0
    aligned_positions = 0
    for query_char, hit_char in zip(query_sequence, hit_sequence):
        if query_char == "-" and hit_char == "-":
            continue
        aligned_positions += 1
        if query_char != "-" and hit_char != "-" and query_char.upper() == hit_char.upper():
            matched += 1

    if aligned_positions == 0:
        aligned_positions = hit.alignment_length
    percent = (100.0 * matched / aligned_positions) if aligned_positions else 0.0
    return f"{matched}/{aligned_positions}({percent:.0f}%)"


def write_legacy_alignment_table(path: Path, hits: list[BlastHit], records_by_id: dict[str, FastaRecord]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write(
            "#Contig_ID\tContig_Seq\tContig_Len\tHit_ID\tHit_Len\tGenus\tDescription\t"
            "Contig_start\tContig_end\tHit_start\tHit_end\tHsp_identity\tE_value\tHsp_strand\n"
        )
        for hit in sorted(hits, key=lambda item: (item.contig_id, item.hit_id, min(item.hit_start, item.hit_end), item.query_start)):
            record = records_by_id.get(hit.contig_id)
            contig_sequence = record.sequence if record is not None else ""
            strand = "-1" if hit.hit_start > hit.hit_end else "1"
            handle.write(
                "\t".join(
                    [
                        hit.contig_id,
                        contig_sequence,
                        str(hit.contig_length),
                        hit.hit_id,
                        str(hit.hit_length),
                        hit.genus,
                        hit.description,
                        str(hit.query_start),
                        str(hit.query_end),
                        str(hit.hit_start),
                        str(hit.hit_end),
                        compute_alignment_identity_text(hit),
                        format(hit.evalue, ".6g"),
                        strand,
                    ]
                )
                + "\n"
            )


def cigar_from_alignment(query_sequence: str, hit_sequence: str) -> str:
    operations: list[tuple[str, int]] = []
    for query_char, hit_char in zip(query_sequence, hit_sequence):
        if query_char == "-" and hit_char == "-":
            continue
        if query_char == "-":
            op = "D"
        elif hit_char == "-":
            op = "I"
        else:
            op = "M"
        if operations and operations[-1][0] == op:
            operations[-1] = (op, operations[-1][1] + 1)
        else:
            operations.append((op, 1))
    return "".join(f"{length}{op}" for op, length in operations) or f"{len(query_sequence.replace('-', ''))}M"


def format_sam_query_sequence(hit: BlastHit) -> str:
    sequence = hit.query_sequence.replace("-", "")
    if hit.hit_start > hit.hit_end:
        return reverse_complement(sequence)
    return sequence


def write_legacy_alignment_sam(path: Path, hits: list[BlastHit]) -> None:
    reference_lengths: dict[str, int] = {}
    for hit in hits:
        reference_lengths.setdefault(hit.hit_id, hit.hit_length)

    with path.open("w", encoding="utf-8", newline="") as handle:
        for hit_id in sorted(reference_lengths):
            handle.write(f"@SQ\tSN:{hit_id}\tLN:{reference_lengths[hit_id]}\n")

        for hit in sorted(hits, key=lambda item: (item.contig_id, item.hit_id, min(item.hit_start, item.hit_end), item.query_start)):
            flag = "16" if hit.hit_start > hit.hit_end else "0"
            query_left_clip = max(0, min(hit.query_start, hit.query_end) - 1)
            query_right_clip = max(0, hit.contig_length - max(hit.query_start, hit.query_end))
            cigar_tokens: list[str] = []
            if query_left_clip:
                cigar_tokens.append(f"{query_left_clip}H")
            cigar_tokens.append(cigar_from_alignment(hit.query_sequence, hit.hit_sequence))
            if query_right_clip:
                cigar_tokens.append(f"{query_right_clip}H")
            if hit.hit_start > hit.hit_end:
                cigar_tokens = list(reversed(cigar_tokens))

            handle.write(
                "\t".join(
                    [
                        hit.contig_id,
                        flag,
                        hit.hit_id,
                        str(min(hit.hit_start, hit.hit_end)),
                        "255",
                        "".join(cigar_tokens),
                        "*",
                        "0",
                        "0",
                        format_sam_query_sequence(hit) or "*",
                        "*",
                        f"AS:i:{int(round(hit.bit_score))}",
                        f"EV:Z:{format(hit.evalue, '.6g')}",
                    ]
                )
                + "\n"
            )


def summarize_hits(
    hits: list[BlastHit],
    contig_depths: dict[str, ContigDepth] | None = None,
    library_size: int = 0,
    coverage_cutoff: float = 0.0,
    depth_cutoff: float = 0.0,
    norm_depth_cutoff: float = 0.0,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    grouped = group_hits_by_reference(hits)
    for key in sorted(grouped):
        analysis, group_id, reference_id, hit_length, genus, description = key
        group_hits = grouped[key]
        best_by_contig = select_best_hits_by_contig(group_hits)

        supporting_contigs = sorted(best_by_contig)
        total_contig_length = sum(best_by_contig[contig_id].contig_length for contig_id in supporting_contigs)
        total_depth = 0
        if contig_depths is not None:
            total_depth = sum(
                contig_depths.get(contig_id, ContigDepth(contig_id, 0, 0, 0, 0.0, 0.0)).total_depth
                for contig_id in supporting_contigs
            )
        covered_bases = merge_covered_bases([(hit.hit_start, hit.hit_end) for hit in group_hits])
        reference_coverage = (covered_bases / hit_length) if hit_length else 0.0
        reference_depth = (total_depth / total_contig_length) if total_contig_length else 0.0
        normalized_depth = (
            (1e6 * total_depth) / (total_contig_length * library_size)
            if total_contig_length and library_size
            else 0.0
        )
        passes_coverage_depth = reference_coverage > coverage_cutoff and (
            reference_depth > depth_cutoff or normalized_depth > norm_depth_cutoff
        )
        rows.append(
            {
                "analysis": analysis,
                "group_id": group_id,
                "reference_id": reference_id,
                "hit_length": hit_length,
                "genus": genus,
                "description": description,
                "contig_count": len(supporting_contigs),
                "supporting_contigs": ",".join(supporting_contigs),
                "total_contig_length": total_contig_length,
                "covered_bases": covered_bases,
                "reference_coverage": round(reference_coverage, 6),
                "reference_depth": round(reference_depth, 6),
                "normalized_depth": round(normalized_depth, 6),
                "mean_percent_identity": round(
                    sum(best_by_contig[contig_id].percent_identity for contig_id in supporting_contigs) / len(supporting_contigs),
                    2,
                ),
                "best_evalue": min(best_by_contig[contig_id].evalue for contig_id in supporting_contigs),
                "max_query_coverage": round(max(best_by_contig[contig_id].query_coverage for contig_id in supporting_contigs), 2),
                "passes_coverage_depth": passes_coverage_depth,
                "passes_redundancy": passes_coverage_depth,
                "passes_filters": passes_coverage_depth,
            }
        )
    return rows


def build_group_support(hits: list[BlastHit]) -> dict[str, dict[str, object]]:
    support: dict[str, dict[str, object]] = {}
    for key, group_hits in group_hits_by_reference(hits).items():
        _, group_id, _, hit_length, _, _ = key
        best_by_contig = select_best_hits_by_contig(group_hits)
        support[group_id] = {
            "hit_length": hit_length,
            "contigs": tuple(sorted(best_by_contig)),
            "ranges": {
                contig_id: tuple(sorted((best_hit.hit_start, best_hit.hit_end)))
                for contig_id, best_hit in best_by_contig.items()
            },
        }
    return support


def prune_redundant_group_ids(
    summary_rows: list[dict[str, object]],
    group_support: dict[str, dict[str, object]],
    diff_ratio: float = 0.25,
    diff_contig_cover: float = 0.5,
    diff_contig_length: int = 100,
) -> tuple[set[str], set[str]]:
    kept_group_ids: list[str] = []
    redundant_group_ids: set[str] = set()
    candidates = sorted(
        (row for row in summary_rows if row["passes_coverage_depth"]),
        key=lambda row: (-int(row["contig_count"]), -int(row["covered_bases"]), str(row["group_id"])),
    )

    for row in candidates:
        group_id = str(row["group_id"])
        support = group_support.get(group_id)
        if support is None:
            kept_group_ids.append(group_id)
            continue

        current_contigs = tuple(str(contig_id) for contig_id in support["contigs"])
        total_contigs = len(current_contigs)
        current_hit_length = int(support["hit_length"])
        current_ranges = support["ranges"]
        is_redundant = False

        for prior_group_id in kept_group_ids:
            prior_support = group_support.get(prior_group_id)
            if prior_support is None:
                continue

            prior_contigs = set(str(contig_id) for contig_id in prior_support["contigs"])
            unique_contigs = [contig_id for contig_id in current_contigs if contig_id not in prior_contigs]
            unique_ratio = (len(unique_contigs) / total_contigs) if total_contigs else 0.0
            unique_coverage_bases = merge_covered_bases([current_ranges[contig_id] for contig_id in unique_contigs])
            unique_coverage_ratio = (unique_coverage_bases / current_hit_length) if current_hit_length else 0.0

            if (
                unique_ratio <= diff_ratio
                and unique_coverage_ratio <= diff_contig_cover
                and unique_coverage_bases <= diff_contig_length
            ):
                redundant_group_ids.add(group_id)
                is_redundant = True
                break

        if not is_redundant:
            kept_group_ids.append(group_id)

    return set(kept_group_ids), redundant_group_ids


def apply_summary_filter_flags(
    summary_rows: list[dict[str, object]],
    kept_group_ids: set[str] | None = None,
) -> list[dict[str, object]]:
    effective_kept_group_ids = kept_group_ids
    if effective_kept_group_ids is None:
        effective_kept_group_ids = {str(row["group_id"]) for row in summary_rows if row["passes_coverage_depth"]}

    updated_rows: list[dict[str, object]] = []
    for row in summary_rows:
        updated_row = dict(row)
        passes_coverage_depth = bool(updated_row["passes_coverage_depth"])
        passes_redundancy = passes_coverage_depth and str(updated_row["group_id"]) in effective_kept_group_ids
        updated_row["passes_redundancy"] = passes_redundancy
        updated_row["passes_filters"] = passes_coverage_depth and passes_redundancy
        updated_rows.append(updated_row)
    return updated_rows


def write_summary_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    field_names = [
        "analysis",
        "group_id",
        "reference_id",
        "hit_length",
        "genus",
        "description",
        "contig_count",
        "supporting_contigs",
        "total_contig_length",
        "covered_bases",
        "reference_coverage",
        "reference_depth",
        "normalized_depth",
        "mean_percent_identity",
        "best_evalue",
        "max_query_coverage",
        "passes_coverage_depth",
        "passes_redundancy",
        "passes_filters",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, field_names, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_contig_depth_tsv(path: Path, contig_depths: dict[str, ContigDepth]) -> None:
    field_names = [
        "contig_id",
        "contig_length",
        "covered_bases",
        "contig_coverage",
        "total_depth",
        "mean_depth",
        "normalized_depth",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, field_names, delimiter="\t")
        writer.writeheader()
        for contig_id in sorted(contig_depths):
            item = contig_depths[contig_id]
            writer.writerow(
                {
                    "contig_id": item.contig_id,
                    "contig_length": item.contig_length,
                    "covered_bases": item.covered_bases,
                    "contig_coverage": round(item.covered_bases / item.contig_length, 6) if item.contig_length else 0.0,
                    "total_depth": item.total_depth,
                    "mean_depth": round(item.mean_depth, 6),
                    "normalized_depth": round(item.normalized_depth, 6),
                }
            )


def filter_records(records: list[FastaRecord], selected_ids: set[str]) -> list[FastaRecord]:
    return [record for record in records if record.seq_id in selected_ids]


def filter_reference_records(records: list[FastaRecord], selected_ids: set[str]) -> list[FastaRecord]:
    return [record for record in records if normalize_seq_id(record.seq_id) in selected_ids]


def select_best_raw_hits(hits: list[BlastHit]) -> dict[str, BlastHit]:
    best_by_contig: dict[str, BlastHit] = {}
    for hit in hits:
        current = best_by_contig.get(hit.contig_id)
        if current is None or hit_sort_key(hit) < hit_sort_key(current):
            best_by_contig[hit.contig_id] = hit
    return best_by_contig


def build_group_identity_stats(hits: list[BlastHit]) -> dict[str, dict[str, float]]:
    grouped: dict[str, list[float]] = defaultdict(list)
    for hit in hits:
        grouped[hit_group_id(hit)].append(hit.percent_identity)

    stats: dict[str, dict[str, float]] = {}
    for group_id, identities in grouped.items():
        stats[group_id] = {
            "mean": sum(identities) / len(identities),
            "max": max(identities),
            "min": min(identities),
        }
    return stats


def format_report_float(value: float, digits: int = 2) -> str:
    return f"{value:.{digits}f}"


def hit_reference_interval(hit: BlastHit) -> tuple[int, int]:
    return (min(hit.hit_start, hit.hit_end), max(hit.hit_start, hit.hit_end))


def describe_hit_strand(hit: BlastHit) -> str:
    return "+" if hit.hit_end >= hit.hit_start else "-"


def identity_color(percent_identity: float) -> str:
    clamped = max(0.0, min(100.0, percent_identity))
    if clamped >= 97.0:
        return "#0f766e"
    if clamped >= 90.0:
        return "#2563eb"
    if clamped >= 75.0:
        return "#d97706"
    return "#b91c1c"


def resolve_ncbi_link(analysis: str, group_id: str, reference_id: str) -> tuple[str, str, str]:
    if analysis == "blastx":
        accession = group_id
        return accession, f"https://www.ncbi.nlm.nih.gov/protein/{html.escape(accession)}", "Protein ID"
    accession = reference_id
    return accession, f"https://www.ncbi.nlm.nih.gov/nuccore/{html.escape(accession)}", "Reference ID"


def sort_hits_by_reference(hits: list[BlastHit]) -> list[BlastHit]:
    return sorted(hits, key=lambda item: (hit_reference_interval(item)[0], hit_reference_interval(item)[1], item.contig_id))


def group_hits_for_reference_panel(hits: list[BlastHit]) -> list[tuple[str, list[BlastHit]]]:
    grouped_hits: dict[str, list[BlastHit]] = defaultdict(list)
    for hit in sort_hits_by_reference(hits):
        grouped_hits[hit.contig_id].append(hit)
    return sorted(
        grouped_hits.items(),
        key=lambda item: (
            min(hit_reference_interval(hit)[0] for hit in item[1]),
            min(hit_reference_interval(hit)[1] for hit in item[1]),
            item[0],
        ),
    )


def describe_panel_track_meta(hits: list[BlastHit]) -> str:
    hsp_count = len(hits)
    strand_values = sorted({describe_hit_strand(hit) for hit in hits})
    strand_summary = strand_values[0] if len(strand_values) == 1 else "mixed"
    identity_values = [hit.percent_identity for hit in hits]
    if min(identity_values) == max(identity_values):
        identity_summary = f"{format_report_float(identity_values[0], 1)}%"
    else:
        identity_summary = f"{format_report_float(min(identity_values), 1)}-{format_report_float(max(identity_values), 1)}%"
    hsp_label = "Alignment" if hsp_count == 1 else "Alignments"
    return f"{hsp_count} {hsp_label} / {strand_summary} / {identity_summary}"


def render_coverage_blocks(reference_length: int, hits: list[BlastHit]) -> str:
    blocks: list[str] = []
    for hit in sort_hits_by_reference(hits):
        ref_start, ref_end = hit_reference_interval(hit)
        left = ((ref_start - 1) / reference_length) * 100
        width = max(((ref_end - ref_start + 1) / reference_length) * 100, 0.35)
        strand_class = "reverse" if describe_hit_strand(hit) == "-" else "forward"
        title = (
            f"{hit.contig_id}: {ref_start}-{ref_end} "
            f"({describe_hit_strand(hit)} strand, {format_report_float(hit.percent_identity, 2)}% identity)"
        )
        blocks.append(
            f"<div class=\"coverage-block {strand_class}\" style=\"left: {left:.3f}%; width: {width:.3f}%;\" "
            f"title=\"{html.escape(title)}\"></div>"
        )
    return "".join(blocks)


def render_reference_track_panel(
    analysis: str,
    group_id: str,
    reference_id: str,
    reference_length: int,
    hits: list[BlastHit],
    anchor_ids: dict[BlastHit, str] | None = None,
) -> str:
    if reference_length <= 0 or not hits:
        return ""

    grouped_tracks = group_hits_for_reference_panel(hits)
    label_width = 140
    track_width = 640
    row_height = 26
    top_padding = 34
    bottom_padding = 18
    svg_width = label_width + track_width + 36
    svg_height = top_padding + bottom_padding + (len(grouped_tracks) * row_height)
    axis_y = 16
    track_left = label_width
    track_right = label_width + track_width
    tick_values = sorted({1, max(reference_length // 4, 1), max(reference_length // 2, 1), max((reference_length * 3) // 4, 1), reference_length})
    link_accession, ncbi_url, link_label = resolve_ncbi_link(analysis, group_id, reference_id)

    parts = [
        "<div class=\"section\" id=\"coverage-map\">",
        "<h2>Reference Coverage Map</h2>",
        "<p class=\"muted\">Reference map with contig tracks, identity-colored segments, and a direct NCBI link for the reference accession.</p>",
        f"<p>{html.escape(link_label)}: <a href=\"{ncbi_url}\" target=\"_blank\" rel=\"noreferrer\"><code>{html.escape(link_accession)}</code></a></p>",
        f"<svg class=\"reference-panel\" viewBox=\"0 0 {svg_width} {svg_height}\" role=\"img\" aria-label=\"Reference coverage panel for {html.escape(group_id)}\">",
        f"<line x1=\"{track_left}\" y1=\"{axis_y}\" x2=\"{track_right}\" y2=\"{axis_y}\" class=\"panel-axis\" />",
    ]

    for tick_value in tick_values:
        x = track_left + ((tick_value - 1) / reference_length) * track_width
        parts.append(f"<line x1=\"{x:.2f}\" y1=\"{axis_y - 5}\" x2=\"{x:.2f}\" y2=\"{axis_y + 5}\" class=\"panel-tick\" />")
        parts.append(f"<text x=\"{x:.2f}\" y=\"10\" text-anchor=\"middle\" class=\"panel-axis-label\">{tick_value}</text>")

    for index, (contig_id, contig_hits) in enumerate(grouped_tracks):
        y_center = top_padding + (index * row_height)
        parts.extend(
            [
                f"<text x=\"{label_width - 8}\" y=\"{y_center + 4}\" text-anchor=\"end\" class=\"panel-label\">{html.escape(contig_id)}</text>",
                f"<line x1=\"{track_left}\" y1=\"{y_center}\" x2=\"{track_right}\" y2=\"{y_center}\" class=\"panel-guide\" />",
            ]
        )
        for hit in contig_hits:
            ref_start, ref_end = hit_reference_interval(hit)
            rect_x = track_left + ((ref_start - 1) / reference_length) * track_width
            rect_width = max(((ref_end - ref_start + 1) / reference_length) * track_width, 3.0)
            color = identity_color(hit.percent_identity)
            anchor = anchor_ids.get(hit) if anchor_ids else None
            title = (
                f"{hit.contig_id}: {ref_start}-{ref_end} "
                f"({describe_hit_strand(hit)} strand, {format_report_float(hit.percent_identity, 2)}% identity, "
                f"E-value {hit.evalue:.6g})"
            )
            if anchor:
                parts.append(f"<a href=\"#{html.escape(anchor)}\">")
            parts.append(
                f"<rect x=\"{rect_x:.2f}\" y=\"{y_center - 7}\" width=\"{rect_width:.2f}\" height=\"14\" rx=\"7\" ry=\"7\" fill=\"{color}\" class=\"panel-segment\">"
                f"<title>{html.escape(title)}</title>"
                "</rect>"
            )
            if anchor:
                parts.append("</a>")
        parts.append(
            f"<text x=\"{track_right + 10}\" y=\"{y_center + 4}\" class=\"panel-meta\">{html.escape(describe_panel_track_meta(contig_hits))}</text>"
        )

    parts.extend(
        [
            "</svg>",
            "<p class=\"muted\">Track color reflects identity: teal `>=97%`, blue `>=90%`, orange `>=75%`, red below that.</p>",
            "</div>",
        ]
    )
    return "".join(parts)


def render_reference_coverage_map(reference_length: int, hits: list[BlastHit]) -> str:
    if reference_length <= 0 or not hits:
        return "<p class=\"muted\">Reference coverage map is unavailable for this reference.</p>"

    return (
        "<div class=\"section\">"
        f"<div class=\"coverage-axis\"><span>1</span><span>{reference_length}</span></div>"
        f"<div class=\"coverage-track\">{render_coverage_blocks(reference_length, hits)}</div>"
        "<p class=\"muted\">Compact overview: each block shows one contig aligned on the reference coordinate axis. Blue indicates forward strand; orange indicates reverse strand.</p>"
        "</div>"
    )


def render_reference_coverage_strip(reference_length: int, hits: list[BlastHit]) -> str:
    if reference_length <= 0 or not hits:
        return "<span class=\"muted\">-</span>"

    forward_count = sum(1 for hit in hits if describe_hit_strand(hit) == "+")
    reverse_count = len(hits) - forward_count
    summary = f"{len(hits)} contigs; +:{forward_count} -:{reverse_count}"
    return (
        f"<div class=\"coverage-track coverage-track-compact\" title=\"{html.escape(summary)}\">"
        f"{render_coverage_blocks(reference_length, hits)}"
        "</div>"
    )


def build_html_document(title: str, body: str) -> str:
    return (
        "<!doctype html>\n"
        "<html lang=\"en\">\n"
        "<head>\n"
        "  <meta charset=\"utf-8\">\n"
        f"  <title>{html.escape(title)}</title>\n"
        "  <style>\n"
        "    body { font-family: Helvetica, Arial, sans-serif; margin: 24px; color: #1f2933; }\n"
        "    h1 { margin-bottom: 8px; }\n"
        "    h2 { margin: 28px 0 8px; font-size: 18px; }\n"
        "    p { line-height: 1.5; }\n"
        "    table { border-collapse: collapse; width: 100%; font-size: 14px; }\n"
        "    th, td { border: 1px solid #d9e2ec; padding: 8px; vertical-align: top; text-align: left; }\n"
        "    th { background: #f0f4f8; }\n"
        "    code { font-family: Menlo, Consolas, monospace; }\n"
        "    .muted { color: #52606d; }\n"
        "    .section { margin-top: 20px; }\n"
        "    .section-nav { margin: 16px 0 20px; font-size: 14px; }\n"
        "    .section-nav a { margin-right: 14px; }\n"
        "    .coverage-axis { display: flex; justify-content: space-between; font-size: 12px; color: #52606d; margin: 6px 0; }\n"
        "    .coverage-track { position: relative; height: 20px; border-radius: 999px; background: #d9e2ec; overflow: hidden; }\n"
        "    .coverage-track-compact { min-width: 180px; height: 14px; }\n"
        "    .coverage-block { position: absolute; top: 0; height: 100%; border-radius: 999px; opacity: 0.9; }\n"
        "    .coverage-block.forward { background: #2f6fed; }\n"
        "    .coverage-block.reverse { background: #d97706; }\n"
        "    .reference-panel { width: 100%; max-width: 860px; height: auto; border: 1px solid #d9e2ec; border-radius: 8px; background: #ffffff; }\n"
        "    .panel-axis, .panel-guide, .panel-tick { stroke: #9fb3c8; stroke-width: 1; }\n"
        "    .panel-guide { stroke-dasharray: 3 3; }\n"
        "    .panel-axis-label, .panel-label, .panel-meta { font-size: 11px; fill: #334e68; font-family: Menlo, Consolas, monospace; }\n"
        "    .panel-segment { opacity: 0.95; }\n"
        "    .alignment-block { margin-top: 16px; padding: 12px; border: 1px solid #d9e2ec; border-radius: 8px; background: #f8fbff; }\n"
        "    .alignment-block table { margin-bottom: 12px; font-size: 13px; }\n"
        "    .alignment-summary { margin: 8px 0; font-size: 13px; line-height: 1.6; }\n"
        "    .alignment-summary-label { font-weight: 600; color: #243b53; margin-right: 8px; }\n"
        "    .alignment-summary-item { display: inline-block; margin-right: 10px; color: #334e68; }\n"
        "    .alignment-segment { margin-top: 14px; padding-top: 14px; border-top: 1px dashed #bcccdc; }\n"
        "    .alignment-segment:first-of-type { margin-top: 0; padding-top: 0; border-top: 0; }\n"
        "    .alignment-header { margin: 8px 0 4px; font-size: 13px; }\n"
        "    .alignment-metrics { margin: 0 0 6px; font-size: 13px; line-height: 1.6; }\n"
        "    .alignment-metric { display: inline-block; margin-right: 10px; }\n"
        "    .alignment-links { margin: 6px 0 10px; font-size: 13px; }\n"
        "    .alignment-text { overflow-x: auto; padding: 12px; border-radius: 6px; background: #102a43; color: #f0f4f8; font-family: Menlo, Consolas, monospace; font-size: 13px; }\n"
        "    .alignment-line { display: block; padding: 2px 8px; border-radius: 2px; white-space: pre; background: #99ffff; color: #102a43; }\n"
        "    .alignment-line + .alignment-line { margin-top: 2px; }\n"
        "    .alignment-spacer { height: 10px; }\n"
        "    .candidate-row td { background: #dcfce7; }\n"
        "    .candidate-flag { display: inline-block; min-width: 44px; padding: 2px 8px; border-radius: 999px; font-size: 12px; font-weight: 600; text-align: center; }\n"
        "    .candidate-yes { background: #15803d; color: #ffffff; }\n"
        "    .candidate-no { background: #d9e2ec; color: #334e68; }\n"
        "  </style>\n"
        "</head>\n"
        "<body>\n"
        f"{body}\n"
        "</body>\n"
        "</html>\n"
    )


def render_html_table(headers: list[str], rows: list[list[str]], row_classes: list[str] | None = None) -> str:
    header_html = "".join(f"<th>{html.escape(header)}</th>" for header in headers)
    body_rows = []
    for index, row in enumerate(rows):
        row_class = ""
        if row_classes is not None and index < len(row_classes) and row_classes[index]:
            row_class = f' class="{html.escape(row_classes[index])}"'
        body_rows.append("<tr" + row_class + ">" + "".join(f"<td>{cell}</td>" for cell in row) + "</tr>")
    return "<table>\n<thead><tr>" + header_html + "</tr></thead>\n<tbody>\n" + "\n".join(body_rows) + "\n</tbody>\n</table>"


def render_undetermined_table(
    headers: list[str],
    rows: list[list[str]],
    row_classes: list[str] | None = None,
    row_ids: list[str] | None = None,
    show_sirna_table: bool = False,
    include_hit_columns: bool = False,
) -> str:
    header_tooltips = {
        "ID": "Contig identifier reported in this analysis.",
        "Contig ID": "Contig identifier reported in this analysis.",
        "Length": "Contig length in bases for this undetermined sequence.",
        "Covered Bases": "Number of contig bases supported by mapped reads.",
        "Contig Coverage": "Percent of the contig span covered by mapped reads.",
        "21-22 (%)": "Percent of 18-33 nt reads assigned to this contig that fall in the 21-22 nt range.",
        "Candidate": "Candidate flag based on the configured 21-22 nt enrichment threshold.",
        "Depth": "Average read depth across the covered contig span.",
        "Mean Depth": "Average read depth across the covered contig span.",
        "Depth (Norm)": "Depth normalized per million input reads for this sample.",
        "Acc#": "Underlying nucleotide reference accession associated with the selected virus-database match.",
        "Match ID": "Primary virus-database match identifier for this undetermined contig.",
        "Best Hit": "Selected virus-database match for this undetermined contig.",
        "Reference ID": "Underlying nucleotide reference accession associated with the selected virus-database match.",
        "Genus": "Reference genus annotation from the database metadata.",
        "BLAST": "BLAST program used for this selected match.",
        "Analysis": "BLAST program used for this selected match.",
        "E value": "Expectation value reported for this selected match.",
        "E-value": "Expectation value reported for this selected match.",
        "Identity": "Percent identity reported for this selected match.",
        "Query Coverage": "Percent of the contig query covered by this selected match.",
        "Description": "Reference description from the database metadata.",
    }
    for length in range(18, 34):
        header_tooltips[str(length)] = (
            f"Number of small-RNA reads of length {length} nt assigned to this undetermined contig."
        )

    body_rows = []
    for index, row in enumerate(rows):
        row_attributes: list[str] = []
        if row_classes is not None and index < len(row_classes) and row_classes[index]:
            row_attributes.append(f'class="{html.escape(row_classes[index])}"')
        if row_ids is not None and index < len(row_ids) and row_ids[index]:
            row_attributes.append(f'id="{html.escape(row_ids[index])}"')
        row_attribute_text = (" " + " ".join(row_attributes)) if row_attributes else ""
        body_rows.append("<tr" + row_attribute_text + ">" + "".join(f"<td>{cell}</td>" for cell in row) + "</tr>")

    has_coverage_columns = "Covered Bases" in headers and "Contig Coverage" in headers
    group_headers = [
        (
            "<th colspan=\"2\" title=\"Contig identifier and assembled length for this undetermined sequence.\">"
            "Contig</th>"
        )
    ]
    if has_coverage_columns:
        group_headers.append(
            (
                "<th colspan=\"2\" title=\"Read-supported coverage metrics for the undetermined contig.\">"
                "Coverage</th>"
            )
        )
    if show_sirna_table:
        group_headers.append(
            (
                "<th colspan=\"17\" title=\"Per-contig read counts across the historical 18-33 nt small-RNA size bins, plus 21-22 nt enrichment.\">"
                "siRNA size distribution</th>"
            )
        )
        group_headers.append(
            "<th colspan=\"1\" title=\"Candidate flag driven by the configured 21-22 nt enrichment threshold.\">Candidate</th>"
        )
        group_headers.append(
            "<th colspan=\"2\" title=\"Read-depth metrics for the assembled undetermined contig.\">Depth</th>"
        )
    else:
        group_headers.append(
            "<th colspan=\"2\" title=\"Read-depth metrics for the assembled undetermined contig.\">Depth</th>"
        )
    if include_hit_columns:
        group_headers.append(
            "<th colspan=\"5\" title=\"Summary columns for the virus-database match on this undetermined contig.\">"
            "Virus-database match</th>"
        )
        group_headers.append(
            "<th colspan=\"3\" title=\"Additional identifiers and alignment metrics for this virus-database match.\">"
            "Match metrics</th>"
        )

    header_row_one = "<tr>" + "".join(group_headers) + "</tr>"
    header_row_two = (
        "<tr>"
        + "".join(
            (
                f"<th title=\"{html.escape(header_tooltips[header])}\">{html.escape(header)}</th>"
                if header in header_tooltips
                else f"<th>{html.escape(header)}</th>"
            )
            for header in headers
        )
        + "</tr>"
    )
    return "<table>\n<thead>" + header_row_one + header_row_two + "</thead>\n<tbody>\n" + "\n".join(body_rows) + "\n</tbody>\n</table>"


def build_summary_row_id(group_id: str) -> str:
    sanitized = "".join(character if character.isalnum() else "-" for character in group_id).strip("-")
    return f"summary-row-{sanitized or 'item'}"


def build_reference_target_map(*analysis_rows: tuple[str, list[dict[str, object]]]) -> dict[tuple[str, str], dict[str, str]]:
    targets: dict[tuple[str, str], dict[str, str]] = {}
    for analysis, rows in analysis_rows:
        for row in rows:
            if not row.get("passes_filters"):
                continue
            group_id = str(row["group_id"])
            targets[(analysis, group_id)] = {
                "summary_href": f"{analysis}.html#{build_summary_row_id(group_id)}",
            }
    return targets


def build_undetermined_row_id(prefix: str, contig_id: str) -> str:
    sanitized = "".join(character if character.isalnum() else "-" for character in contig_id).strip("-")
    return f"undetermined-{prefix}-row-{sanitized or 'item'}"


def build_undetermined_contig_cell(contig_id: str, row_id: str | None = None) -> str:
    contig_label = f"<code>{html.escape(contig_id)}</code>"
    if row_id is None:
        return contig_label
    return f"<a href=\"#{html.escape(row_id)}\">{contig_label}</a>"


def render_reference_summary_table(
    headers: list[str],
    rows: list[list[str]],
    row_ids: list[str] | None = None,
) -> str:
    header_tooltips = {
        "Reference": "Reference accession or grouped identifier for this report row.",
        "Reference ID": "Underlying nucleotide reference accession associated with this reference group.",
        "Length": "Reference sequence length in bases or amino acids for this group.",
        "Coverage": "Covered bases and percentage of the reference span supported by grouped contigs.",
        "Map": "Compact coordinate strip showing where grouped contigs align on the reference axis.",
        "#Contig": "Number of contigs assigned to this reference group.",
        "Depth": "Average depth across the covered reference span.",
        "Depth (Norm)": "Depth normalized per million input reads for this sample.",
        "%Identity": "Mean percent identity across contig alignments for this reference group.",
        "%Iden Max": "Maximum percent identity among contig alignments for this reference group.",
        "%Iden Min": "Minimum percent identity among contig alignments for this reference group.",
        "Genus": "Reference genus annotation from the database metadata.",
        "Description": "Reference description from the database metadata.",
    }
    header_row_one = (
        "<tr>"
        "<th colspan=\"2\">Reference</th>"
        "<th colspan=\"6\">Coverage</th>"
        "<th colspan=\"3\">Identity</th>"
        "<th colspan=\"2\">Annotation</th>"
        "</tr>"
    )
    header_row_two = (
        "<tr>"
        + "".join(
            (
                f"<th title=\"{html.escape(header_tooltips[header])}\">{html.escape(header)}</th>"
                if header in header_tooltips
                else f"<th>{html.escape(header)}</th>"
            )
            for header in headers
        )
        + "</tr>"
    )
    body_rows = []
    for index, row in enumerate(rows):
        row_id = ""
        if row_ids is not None and index < len(row_ids):
            row_id = f' id="{html.escape(row_ids[index])}"'
        body_rows.append("<tr" + row_id + ">" + "".join(f"<td>{cell}</td>" for cell in row) + "</tr>")
    return "<table>\n<thead>" + header_row_one + header_row_two + "</thead>\n<tbody>\n" + "\n".join(body_rows) + "\n</tbody>\n</table>"


def render_reference_detail_hits_table(
    headers: list[str],
    rows: list[list[str]],
    row_ids: list[str] | None = None,
) -> str:
    header_tooltips = {
        "Order": "Stable order of alignments on this reference detail page.",
        "Contig ID": "Contig identifier reported in this analysis.",
        "Contig Length": "Contig length in bases or amino acids for this alignment.",
        "Query Start": "Start coordinate of the aligned segment on the contig sequence.",
        "Query End": "End coordinate of the aligned segment on the contig sequence.",
        "Reference Start": "Start coordinate of the aligned segment on the reference sequence.",
        "Reference End": "End coordinate of the aligned segment on the reference sequence.",
        "Strand": "Reference orientation for this alignment.",
        "Identity": "Percent identity reported for this alignment.",
        "E-value": "Expectation value reported for this alignment.",
        "Query Coverage": "Percent of the contig query covered by this alignment.",
        "Description": "Reference description from the database metadata.",
    }
    header_row_one = (
        "<tr>"
        "<th colspan=\"3\">Contig</th>"
        "<th colspan=\"5\">Alignment Coordinates</th>"
        "<th colspan=\"3\">Metrics</th>"
        "<th colspan=\"1\">Annotation</th>"
        "</tr>"
    )
    header_row_two = (
        "<tr>"
        + "".join(
            (
                f"<th title=\"{html.escape(header_tooltips[header])}\">{html.escape(header)}</th>"
                if header in header_tooltips
                else f"<th>{html.escape(header)}</th>"
            )
            for header in headers
        )
        + "</tr>"
    )
    body_rows = []
    for index, row in enumerate(rows):
        row_id = ""
        if row_ids is not None and index < len(row_ids):
            row_id = f' id="{html.escape(row_ids[index])}"'
        body_rows.append("<tr" + row_id + ">" + "".join(f"<td>{cell}</td>" for cell in row) + "</tr>")
    return "<table>\n<thead>" + header_row_one + header_row_two + "</thead>\n<tbody>\n" + "\n".join(body_rows) + "\n</tbody>\n</table>"


def render_alignment_block_hits_table(rows: list[list[str]], row_ids: list[str]) -> str:
    headers = [
        "Order",
        "Query ID",
        "Query Start",
        "Query End",
        "Reference Start",
        "Reference End",
        "Identity",
        "E-value",
        "Alignment Length",
        "Bit Score",
        "Query Coverage",
        "Strand",
    ]
    header_tooltips = {
        "Order": "Stable order of alignments within this contig block.",
        "Query ID": "Contig identifier for this alignment row.",
        "Query Start": "Start coordinate of the aligned segment on the contig sequence.",
        "Query End": "End coordinate of the aligned segment on the contig sequence.",
        "Reference Start": "Start coordinate of the aligned segment on the reference sequence.",
        "Reference End": "End coordinate of the aligned segment on the reference sequence.",
        "Identity": "Percent identity reported for this alignment.",
        "E-value": "Expectation value reported for this alignment.",
        "Alignment Length": "Aligned length for this contig/reference segment.",
        "Bit Score": "Bit score for this alignment.",
        "Query Coverage": "Percent of the contig query covered by this alignment.",
        "Strand": "Reference orientation for this alignment.",
    }
    header_row_one = (
        "<tr>"
        "<th colspan=\"1\">Navigation</th>"
        "<th colspan=\"3\">Query</th>"
        "<th colspan=\"2\">Reference</th>"
        "<th colspan=\"6\">Metrics</th>"
        "</tr>"
    )
    header_row_two = (
        "<tr>"
        + "".join(
            f"<th title=\"{html.escape(header_tooltips[header])}\">{html.escape(header)}</th>"
            for header in headers
        )
        + "</tr>"
    )
    body_rows = []
    for index, row in enumerate(rows):
        row_id = f' id="{html.escape(row_ids[index])}"' if index < len(row_ids) else ""
        body_rows.append("<tr" + row_id + ">" + "".join(f"<td>{cell}</td>" for cell in row) + "</tr>")
    return "<table>\n<thead>" + header_row_one + header_row_two + "</thead>\n<tbody>\n" + "\n".join(body_rows) + "\n</tbody>\n</table>"


def render_alignment_block_index_table(rows: list[list[str]], row_ids: list[str] | None = None) -> str:
    headers = ["Contig", "Alignments", "Order", "Reference Span", "Identity Range", "Links"]
    header_tooltips = {
        "Contig": "Contig identifier for this contig block.",
        "Alignments": "Number of alignments grouped into this contig block.",
        "Order": "Stable order range covered by this contig block.",
        "Reference Span": "Reference span covered by this contig block.",
        "Identity Range": "Percent-identity range across the alignments in this contig block.",
        "Links": "Quick jumps to this contig block, its first or last alignment, and the matching contig row.",
    }
    header_row_one = (
        "<tr>"
        "<th colspan=\"1\" title=\"Anchor for each contig block listed in this index.\">Contig Block</th>"
        "<th colspan=\"3\" title=\"Alignment count, stable order range, and covered reference span for this contig block.\">Block Summary</th>"
        "<th colspan=\"1\" title=\"Percent-identity range across the alignments in this contig block.\">Identity</th>"
        "<th colspan=\"1\" title=\"Direct jumps to this contig block, its alignments, and the matching contig row.\">Navigation</th>"
        "</tr>"
    )
    header_row_two = (
        "<tr>"
        + "".join(
            f"<th title=\"{html.escape(header_tooltips[header])}\">{html.escape(header)}</th>"
            for header in headers
        )
        + "</tr>"
    )
    body_rows = []
    for index, row in enumerate(rows):
        row_id = ""
        if row_ids is not None and index < len(row_ids):
            row_id = f' id="{html.escape(row_ids[index])}"'
        body_rows.append("<tr" + row_id + ">" + "".join(f"<td>{cell}</td>" for cell in row) + "</tr>")
    return "<table>\n<thead>" + header_row_one + header_row_two + "</thead>\n<tbody>\n" + "\n".join(body_rows) + "\n</tbody>\n</table>"


def build_alignment_midline(query_sequence: str, hit_sequence: str) -> str:
    midline: list[str] = []
    for query_char, hit_char in zip(query_sequence, hit_sequence):
        if query_char == hit_char and query_char != "-":
            midline.append("|")
        elif query_char == "-" or hit_char == "-":
            midline.append(" ")
        else:
            midline.append(".")
    return "".join(midline)


def split_alignment_segments(text: str, width: int = ALIGNMENT_WRAP_WIDTH) -> list[str]:
    if not text:
        return []
    return [text[index : index + width] for index in range(0, len(text), width)]


def consume_alignment_coordinates(segment: str, current: int, step: int) -> tuple[int, int, int]:
    non_gap_count = sum(1 for char in segment if char != "-")
    if non_gap_count <= 0:
        return current, current, current
    end = current + (step * (non_gap_count - 1))
    next_value = end + step
    return current, end, next_value


def format_alignment_line(label: str, start: int, segment: str, end: int, coordinate_width: int) -> str:
    return f"{label:<5} {start:>{coordinate_width}} {segment} {end:>{coordinate_width}}"


def format_alignment_midline(segment: str, coordinate_width: int) -> str:
    return f"{'Match':<5} {'':>{coordinate_width}} {segment}"


def render_alignment_text(hit: BlastHit) -> str:
    query_sequence = hit.query_sequence.strip()
    hit_sequence = hit.hit_sequence.strip()
    if not query_sequence or not hit_sequence:
        return ""

    midline = build_alignment_midline(query_sequence, hit_sequence)
    query_segments = split_alignment_segments(query_sequence)
    midline_segments = split_alignment_segments(midline)
    hit_segments = split_alignment_segments(hit_sequence)
    query_step = 1 if hit.query_end >= hit.query_start else -1
    hit_step = 1 if hit.hit_end >= hit.hit_start else -1
    coordinate_width = max(
        len(str(value))
        for value in (
            hit.query_start,
            hit.query_end,
            hit.hit_start,
            hit.hit_end,
        )
    )
    query_position = hit.query_start
    hit_position = hit.hit_start
    lines: list[str] = []
    for index, query_segment in enumerate(query_segments):
        if index:
            lines.append("")
        hit_segment = hit_segments[index]
        query_start, query_end, query_position = consume_alignment_coordinates(query_segment, query_position, query_step)
        hit_start, hit_end, hit_position = consume_alignment_coordinates(hit_segment, hit_position, hit_step)
        lines.append(format_alignment_line("Query", query_start, query_segment, query_end, coordinate_width))
        lines.append(format_alignment_midline(midline_segments[index], coordinate_width))
        lines.append(format_alignment_line("Hit", hit_start, hit_segment, hit_end, coordinate_width))
    return "\n".join(lines)


def render_alignment_text_html(hit: BlastHit) -> str:
    alignment_text = render_alignment_text(hit)
    if not alignment_text:
        return ""

    html_lines: list[str] = []
    for line in alignment_text.splitlines():
        if not line:
            html_lines.append("<div class=\"alignment-spacer\"></div>")
            continue
        label = line.split(None, 1)[0].lower()
        html_lines.append(
            "<div "
            f"class=\"alignment-line alignment-line-{html.escape(label)}\">"
            f"{html.escape(line)}"
            "</div>"
        )
    return "<div class=\"alignment-text\">" + "".join(html_lines) + "</div>"


def build_alignment_anchor_ids(hits: list[BlastHit]) -> tuple[dict[BlastHit, str], dict[str, str]]:
    hit_anchor_ids: dict[BlastHit, str] = {}
    contig_anchor_ids: dict[str, str] = {}
    next_hit_index = 1
    next_contig_index = 1
    for hit in hits:
        hit_anchor_ids[hit] = f"alignment-hit-{next_hit_index}"
        next_hit_index += 1
        if hit.contig_id not in contig_anchor_ids:
            contig_anchor_ids[hit.contig_id] = f"alignment-group-{next_contig_index}"
            next_contig_index += 1
    return hit_anchor_ids, contig_anchor_ids


def order_alignment_block_contigs(hits: list[BlastHit]) -> tuple[dict[str, list[BlastHit]], list[str]]:
    grouped_hits: dict[str, list[BlastHit]] = defaultdict(list)
    for hit in hits:
        grouped_hits[hit.contig_id].append(hit)
    ordered_contigs = sorted(
        grouped_hits,
        key=lambda contig_id: (
            min(hit_reference_interval(hit)[0] for hit in grouped_hits[contig_id]),
            min(hit_reference_interval(hit)[1] for hit in grouped_hits[contig_id]),
            contig_id,
        ),
    )
    return grouped_hits, ordered_contigs


def build_alignment_index_row_ids(hits: list[BlastHit]) -> dict[str, str]:
    _, ordered_contigs = order_alignment_block_contigs(hits)
    return {
        contig_id: f"alignment-index-row-{index}"
        for index, contig_id in enumerate(ordered_contigs, start=1)
    }


def render_alignment_block_index(
    hits: list[BlastHit],
    contig_anchor_ids: dict[str, str],
    hit_anchor_ids: dict[BlastHit, str],
    detail_hit_row_ids: dict[BlastHit, str],
    index_row_ids: dict[str, str],
) -> str:
    grouped_hits, ordered_contigs = order_alignment_block_contigs(hits)
    if not ordered_contigs:
        return ""

    rows: list[list[str]] = []
    order_start = 1
    for contig_id in ordered_contigs:
        contig_hits = grouped_hits[contig_id]
        reference_starts = [hit_reference_interval(hit)[0] for hit in contig_hits]
        reference_ends = [hit_reference_interval(hit)[1] for hit in contig_hits]
        identity_values = [hit.percent_identity for hit in contig_hits]
        order_end = order_start + len(contig_hits) - 1
        first_hit = contig_hits[0]
        last_hit = contig_hits[-1]
        order_label = (
            f"<a href=\"#{html.escape(detail_hit_row_ids[contig_hits[0]])}\">{html.escape(str(order_start))}</a>"
            if order_start == order_end
            else (
                f"<a href=\"#{html.escape(detail_hit_row_ids[contig_hits[0]])}\">{html.escape(str(order_start))}</a>-"
                f"<a href=\"#{html.escape(detail_hit_row_ids[contig_hits[-1]])}\">{html.escape(str(order_end))}</a>"
            )
        )
        rows.append(
            [
                f"<a href=\"#{html.escape(contig_anchor_ids[contig_id])}\"><code>{html.escape(contig_id)}</code></a>",
                (
                    f"<a href=\"#{html.escape(contig_anchor_ids[contig_id])}\">"
                    f"{html.escape(str(len(contig_hits)))}"
                    "</a>"
                ),
                order_label,
                (
                    f"<a href=\"#{html.escape(contig_anchor_ids[contig_id])}\">"
                    f"{html.escape(str(min(reference_starts)))}-{html.escape(str(max(reference_ends)))}"
                    "</a>"
                ),
                (
                    (
                        f"<a href=\"#{html.escape(hit_anchor_ids[first_hit])}-text\">"
                        f"{html.escape(format_report_float(identity_values[0], 2))}%"
                        "</a>"
                    )
                    if min(identity_values) == max(identity_values)
                    else (
                        f"<a href=\"#{html.escape(hit_anchor_ids[first_hit])}-text\">"
                        f"{html.escape(format_report_float(min(identity_values), 2))}-"
                        f"{html.escape(format_report_float(max(identity_values), 2))}%"
                        "</a>"
                    )
                ),
                (
                    (
                        f"<a href=\"#{html.escape(contig_anchor_ids[contig_id])}\">Contig Block</a> | "
                        f"<a href=\"#{html.escape(hit_anchor_ids[first_hit])}-text\">Alignment</a> | "
                        f"<a href=\"#{html.escape(detail_hit_row_ids[first_hit])}\">Contig Row</a>"
                    )
                    if len(contig_hits) == 1
                    else (
                        f"<a href=\"#{html.escape(contig_anchor_ids[contig_id])}\">Contig Block</a> | "
                        f"<a href=\"#{html.escape(hit_anchor_ids[first_hit])}-text\">First Alignment</a> | "
                        f"<a href=\"#{html.escape(hit_anchor_ids[last_hit])}-text\">Last Alignment</a> | "
                        f"<a href=\"#{html.escape(detail_hit_row_ids[first_hit])}\">Contig Row</a>"
                    )
                ),
            ]
        )
        order_start = order_end + 1

    return render_alignment_block_index_table(rows, [index_row_ids[contig_id] for contig_id in ordered_contigs])


def render_alignment_block(
    contig_id: str,
    hits: list[BlastHit],
    order_start: int,
    block_anchor_id: str,
    index_row_anchor_id: str,
    hit_anchor_ids: dict[BlastHit, str],
    detail_hit_row_ids: dict[BlastHit, str],
    previous_block_anchor_id: str | None = None,
    next_block_anchor_id: str | None = None,
) -> str:
    header_rows: list[list[str]] = []
    header_row_ids: list[str] = []
    alignment_sections: list[str] = []
    reference_starts: list[int] = []
    reference_ends: list[int] = []
    identity_values: list[float] = []
    query_coverage_values: list[float] = []
    contig_length = hits[0].contig_length if hits else 0
    for offset, hit in enumerate(hits):
        alignment_html = render_alignment_text_html(hit)
        if not alignment_html:
            continue
        ref_start, ref_end = hit_reference_interval(hit)
        reference_starts.append(ref_start)
        reference_ends.append(ref_end)
        identity_values.append(hit.percent_identity)
        query_coverage_values.append(hit.query_coverage)
        detail_row_anchor = detail_hit_row_ids[hit]
        text_anchor = f"{hit_anchor_ids[hit]}-text"
        header_rows.append(
            [
                f"<a href=\"#{html.escape(detail_row_anchor)}\">{html.escape(str(order_start + offset))}</a>",
                f"<a href=\"#{html.escape(text_anchor)}\"><code>{html.escape(hit.contig_id)}</code></a>",
                f"<a href=\"#{html.escape(text_anchor)}\">{html.escape(str(hit.query_start))}</a>",
                f"<a href=\"#{html.escape(text_anchor)}\">{html.escape(str(hit.query_end))}</a>",
                f"<a href=\"#{html.escape(text_anchor)}\">{html.escape(str(ref_start))}</a>",
                f"<a href=\"#{html.escape(text_anchor)}\">{html.escape(str(ref_end))}</a>",
                f"<a href=\"#{html.escape(text_anchor)}\">{html.escape(format_report_float(hit.percent_identity, 2))}</a>",
                f"<a href=\"#{html.escape(text_anchor)}\">{html.escape(f'{hit.evalue:.6g}')}</a>",
                f"<a href=\"#{html.escape(text_anchor)}\">{html.escape(str(hit.alignment_length))}</a>",
                f"<a href=\"#{html.escape(text_anchor)}\">{html.escape(format_report_float(hit.bit_score, 1))}</a>",
                f"<a href=\"#{html.escape(text_anchor)}\">{html.escape(format_report_float(hit.query_coverage, 2))}</a>",
                f"<a href=\"#{html.escape(text_anchor)}\">{html.escape(describe_hit_strand(hit))}</a>",
            ]
        )
        header_row_ids.append(hit_anchor_ids[hit])
        segment_title_html = ""
        if len(hits) > 1:
            segment_title_html = (
                "<p>"
                f"<strong>Alignment Segment {html.escape(str(offset + 1))}</strong> "
                f"<span class=\"muted\">(Order <a href=\"#{html.escape(detail_row_anchor)}\">{html.escape(str(order_start + offset))}</a>)</span>"
                "</p>"
            )
        alignment_header_html = (
            "<p class=\"alignment-header\">"
            "<strong>Alignment:</strong> "
            f"<span class=\"muted\">Query {html.escape(str(hit.query_start))}-{html.escape(str(hit.query_end))} | "
            f"Reference {html.escape(str(ref_start))}-{html.escape(str(ref_end))} | "
            f"Strand {html.escape(describe_hit_strand(hit))}</span>"
            "</p>"
        )
        alignment_metrics_html = (
            "<p class=\"alignment-metrics muted\">"
            "<strong>Metrics:</strong> "
            f"<span class=\"alignment-metric\"><strong>Identity:</strong> {html.escape(format_report_float(hit.percent_identity, 2))}%</span>"
            f"<span class=\"alignment-metric\"><strong>E-value:</strong> {html.escape(f'{hit.evalue:.6g}')}</span>"
            f"<span class=\"alignment-metric\"><strong>Bit Score:</strong> {html.escape(format_report_float(hit.bit_score, 1))}</span>"
            f"<span class=\"alignment-metric\"><strong>Alignment Length:</strong> {html.escape(str(hit.alignment_length))}</span>"
            "</p>"
        )
        sibling_links: list[str] = []
        if len(hits) > 1:
            if offset > 0:
                previous_anchor = f"{hit_anchor_ids[hits[offset - 1]]}-text"
                sibling_links.append(f"<a href=\"#{html.escape(previous_anchor)}\">Previous Segment</a>")
            if offset < len(hits) - 1:
                next_anchor = f"{hit_anchor_ids[hits[offset + 1]]}-text"
                sibling_links.append(f"<a href=\"#{html.escape(next_anchor)}\">Next Segment</a>")
        alignment_sections.append(
            (
                f"<div class=\"alignment-segment\" id=\"{html.escape(hit_anchor_ids[hit])}-text\">"
                + segment_title_html
                + alignment_header_html
                + alignment_metrics_html
                + "<p class=\"alignment-links muted\">"
                + f"<a href=\"#{html.escape(detail_row_anchor)}\">Contig Hit Row</a> | "
                + f"<a href=\"#{html.escape(hit_anchor_ids[hit])}\">Block Row</a> | "
                + (" | ".join(sibling_links) + " | " if sibling_links else "")
                + f"<a href=\"#{html.escape(index_row_anchor_id)}\">Index Entry</a> | "
                + "<a href=\"#alignment-index\">Alignment Index</a> | "
                + "<a href=\"#coverage-map\">Coverage Map</a> | "
                + "<a href=\"#detail-top\">Top</a>"
                + "</p>"
                + f"{alignment_html}"
                + "</div>"
            )
        )

    if not header_rows or not alignment_sections:
        return ""

    first_alignment_text_anchor = f"{header_row_ids[0]}-text"
    order_end = order_start + len(header_rows) - 1
    order_range_value = (
        f"<a href=\"#{html.escape(detail_hit_row_ids[hits[0]])}\">{html.escape(str(order_start))}</a>"
        if order_start == order_end
        else (
            f"<a href=\"#{html.escape(detail_hit_row_ids[hits[0]])}\">{html.escape(str(order_start))}</a>-"
            f"<a href=\"#{html.escape(detail_hit_row_ids[hits[-1]])}\">{html.escape(str(order_end))}</a>"
        )
    )
    metadata_items = [
        ("Contig Length", html.escape(str(contig_length))),
        ("Alignments", html.escape(str(len(header_rows)))),
        ("Order Range", order_range_value),
        (
            "Reference Span",
            f"{html.escape(str(min(reference_starts)))}-{html.escape(str(max(reference_ends)))}",
        ),
        (
            "Identity Range",
            f"{html.escape(format_report_float(min(identity_values), 2))}-"
            f"{html.escape(format_report_float(max(identity_values), 2))}%",
        ),
        (
            "Query Coverage Range",
            f"{html.escape(format_report_float(min(query_coverage_values), 2))}-"
            f"{html.escape(format_report_float(max(query_coverage_values), 2))}%",
        ),
    ]
    metadata_html = (
        "<p class=\"alignment-summary\">"
        "<span class=\"alignment-summary-label\">Summary:</span>"
        + "".join(
            (
                "<span class=\"alignment-summary-item\">"
                f"<strong>{label}:</strong> {value}"
                "</span>"
            )
            for label, value in metadata_items
        )
        + "</p>"
    )

    block_links = [
        "<a href=\"#coverage-map\">Coverage Map</a>",
        "<a href=\"#hit-table\">Contig Hits</a>",
        f"<a href=\"#{html.escape(index_row_anchor_id)}\">Index Entry</a>",
        "<a href=\"#alignment-index\">Alignment Index</a>",
    ]
    if previous_block_anchor_id is not None:
        block_links.append(f"<a href=\"#{html.escape(previous_block_anchor_id)}\">Previous Contig Block</a>")
    if next_block_anchor_id is not None:
        block_links.append(f"<a href=\"#{html.escape(next_block_anchor_id)}\">Next Contig Block</a>")
    block_links.append("<a href=\"#detail-top\">Top</a>")

    return (
        f"<div class=\"alignment-block\" id=\"{html.escape(block_anchor_id)}\">"
        f"<p><strong><a href=\"#{html.escape(first_alignment_text_anchor)}\">{html.escape(contig_id)}</a></strong> "
        "<span class=\"muted\">("
        + " | ".join(block_links)
        + ")</span></p>"
        f"{metadata_html}"
        "<p class=\"muted\">Retained alignments for this contig.</p>"
        + render_alignment_block_hits_table(header_rows, header_row_ids)
        + "".join(alignment_sections)
        + "</div>"
    )


def render_alignment_blocks(
    hits: list[BlastHit],
    hit_anchor_ids: dict[BlastHit, str],
    contig_anchor_ids: dict[str, str],
    detail_hit_row_ids: dict[BlastHit, str],
    index_row_ids: dict[str, str],
) -> str:
    grouped_hits, ordered_contigs = order_alignment_block_contigs(hits)
    blocks: list[str] = []
    order_start = 1
    for contig_index, contig_id in enumerate(ordered_contigs):
        contig_hits = grouped_hits[contig_id]
        previous_block_anchor_id = (
            contig_anchor_ids[ordered_contigs[contig_index - 1]]
            if contig_index > 0
            else None
        )
        next_block_anchor_id = (
            contig_anchor_ids[ordered_contigs[contig_index + 1]]
            if contig_index < len(ordered_contigs) - 1
            else None
        )
        block = render_alignment_block(
            contig_id,
            contig_hits,
            order_start,
            contig_anchor_ids[contig_id],
            index_row_ids[contig_id],
            hit_anchor_ids,
            detail_hit_row_ids,
            previous_block_anchor_id=previous_block_anchor_id,
            next_block_anchor_id=next_block_anchor_id,
        )
        if block:
            blocks.append(block)
            order_start += len(contig_hits)
    blocks = [block for block in blocks if block]
    if not blocks:
        return ""

    return (
        "<div class=\"section\" id=\"alignment-blocks\">"
        "<h2>Detailed Alignments</h2>"
        "<p class=\"muted\">Alignment text for each contig segment is shown below.</p>"
        + "".join(blocks)
        + "</div>"
    )


def build_reference_report_html(
    analysis: str,
    summary_rows: list[dict[str, object]],
    final_hits: list[BlastHit],
    reference_fasta_name: str,
    detail_dir_name: str | None = None,
) -> str:
    title = f"{analysis.upper()} Reference Summary"
    kept_rows = sorted(
        (row for row in summary_rows if row["passes_filters"]),
        key=lambda row: (-int(row["contig_count"]), -int(row["covered_bases"]), str(row["group_id"])),
    )
    identity_stats = build_group_identity_stats(final_hits)
    grouped_hits: dict[str, list[BlastHit]] = defaultdict(list)
    for hit in final_hits:
        grouped_hits[hit_group_id(hit)].append(hit)
    body_parts = [
        f"<h1 id=\"detail-top\">{html.escape(title)}</h1>",
        "<p class=\"muted\">Summary report with grouped coverage and identity metrics for retained references.</p>",
        f"<p>Reference FASTA: <code>{html.escape(reference_fasta_name)}</code></p>",
        "<p class=\"section-nav\"><a href=\"#hit-table\">Reference Table</a> | <a href=\"#detail-top\">Top</a></p>",
    ]

    if not kept_rows:
        body_parts.append("<p>No references passed the final filters.</p>")
        return build_html_document(title, "\n".join(body_parts))

    headers = [
        "Reference",
        "Reference ID",
        "Length",
        "Coverage",
        "Map",
        "#Contig",
        "Depth",
        "Depth (Norm)",
        "%Identity",
        "%Iden Max",
        "%Iden Min",
        "Genus",
        "Description",
    ]
    rows: list[list[str]] = []
    row_ids: list[str] = []
    for row in kept_rows:
        group_id = str(row["group_id"])
        row_ids.append(build_summary_row_id(group_id))
        stats = identity_stats.get(group_id, {})
        group_hits = grouped_hits.get(group_id, [])
        detail_href = f"{html.escape(detail_dir_name)}/{html.escape(group_id)}.html" if detail_dir_name else ""
        coverage_text = html.escape(f"{row['covered_bases']} ({float(row['reference_coverage']) * 100:.1f}%)")
        contig_count_text = html.escape(str(row["contig_count"]))
        description_text = html.escape(str(row["description"]))
        reference_id_text = html.escape(str(row["reference_id"]))
        group_id_text = html.escape(group_id)
        rows.append(
            [
                (
                    f"<a href=\"{detail_href}\"><code>{group_id_text}</code></a>"
                    if detail_dir_name
                    else f"<code>{group_id_text}</code>"
                ),
                (
                    f"<a href=\"{detail_href}\"><code>{reference_id_text}</code></a>"
                    if detail_dir_name
                    else f"<code>{reference_id_text}</code>"
                ),
                html.escape(str(row["hit_length"])),
                (
                    f"<a href=\"{detail_href}\">{coverage_text}</a>"
                    if detail_dir_name
                    else coverage_text
                ),
                (
                    f"<a href=\"{detail_href}\">{render_reference_coverage_strip(int(row['hit_length']), group_hits)}</a>"
                    if detail_dir_name
                    else render_reference_coverage_strip(int(row["hit_length"]), group_hits)
                ),
                (
                    f"<a href=\"{detail_href}\">{contig_count_text}</a>"
                    if detail_dir_name
                    else contig_count_text
                ),
                html.escape(format_report_float(float(row["reference_depth"]), 1)),
                html.escape(format_report_float(float(row["normalized_depth"]), 1)),
                html.escape(format_report_float(float(stats.get("mean", row["mean_percent_identity"])), 2)),
                html.escape(format_report_float(float(stats.get("max", row["mean_percent_identity"])), 2)),
                html.escape(format_report_float(float(stats.get("min", row["mean_percent_identity"])), 2)),
                html.escape(str(row["genus"])),
                (
                    f"<a href=\"{detail_href}\">{description_text}</a>"
                    if detail_dir_name
                    else description_text
                ),
            ]
        )

    body_parts.append("<div class=\"section\" id=\"hit-table\"><h2>Reference Table</h2></div>")
    body_parts.append(render_reference_summary_table(headers, rows, row_ids=row_ids))
    return build_html_document(title, "\n".join(body_parts))


def build_reference_detail_html(
    analysis: str,
    group_id: str,
    reference_id: str,
    hits: list[BlastHit],
    summary_row: dict[str, object] | None = None,
    summary_page_name: str | None = None,
) -> str:
    title = f"{analysis.upper()} Reference Detail: {group_id}"
    ordered_hits = sorted(hits, key=lambda item: (hit_reference_interval(item)[0], hit_reference_interval(item)[1], item.contig_id))
    hit_anchor_ids, contig_anchor_ids = build_alignment_anchor_ids(ordered_hits)
    index_row_ids = build_alignment_index_row_ids(ordered_hits)
    detail_hit_row_ids = {hit: f"hit-row-{index}" for index, hit in enumerate(ordered_hits, start=1)}
    body_parts = [
        f"<h1 id=\"detail-top\">{html.escape(title)}</h1>",
        f"<p class=\"muted\">Retained contig alignments for <code>{html.escape(group_id)}</code>.</p>",
    ]

    if summary_row is not None:
        reference_length = int(summary_row.get("hit_length", 0))
        detail_lines = []
        if analysis == "blastx":
            detail_lines.append(f"Protein ID: <code>{html.escape(group_id)}</code>")
            detail_lines.append(f"Reference ID: <code>{html.escape(reference_id)}</code>")
        else:
            detail_lines.append(f"Reference ID: <code>{html.escape(reference_id)}</code>")
        detail_lines.extend(
            [
                f"Reference Length: {html.escape(str(reference_length))}",
                f"Coverage: {html.escape(str(summary_row['covered_bases']))} bases ({float(summary_row['reference_coverage']) * 100:.1f}%)",
                f"Contigs: {html.escape(str(summary_row['contig_count']))}",
                f"Depth: {html.escape(format_report_float(float(summary_row['reference_depth']), 1))}",
                f"Depth (Norm): {html.escape(format_report_float(float(summary_row['normalized_depth']), 1))}",
            ]
        )
        body_parts.append("<p>" + "<br>".join(detail_lines) + "</p>")
    if summary_page_name:
        summary_anchor = ""
        if summary_row is not None:
            summary_anchor = "#" + build_summary_row_id(str(summary_row.get("group_id", group_id)))
        body_parts.append(
            f"<p><a href=\"../{html.escape(summary_page_name)}{summary_anchor}\">Back to {html.escape(analysis.upper())} summary</a></p>"
        )
        body_parts.append(
            "<p>"
            "<a href=\"#coverage-map\">Coverage Map</a> | "
            "<a href=\"#hit-table\">Contig Hits</a> | "
            "<a href=\"#alignment-index\">Alignment Index</a> | "
            "<a href=\"#alignment-blocks\">Detailed Alignments</a>"
            "</p>"
        )

    if not hits:
        body_parts.append("<p>No contig alignments were recorded for this reference.</p>")
        return build_html_document(title, "\n".join(body_parts))

    if summary_row is not None:
        reference_length = int(summary_row.get("hit_length", 0))
        body_parts.append(render_reference_track_panel(analysis, group_id, reference_id, reference_length, ordered_hits, hit_anchor_ids))
        body_parts.append(render_reference_coverage_map(reference_length, ordered_hits))

    headers = [
        "Order",
        "Contig ID",
        "Contig Length",
        "Query Start",
        "Query End",
        "Reference Start",
        "Reference End",
        "Strand",
        "Identity",
        "E-value",
        "Query Coverage",
        "Description",
    ]
    rows: list[list[str]] = []
    for index, hit in enumerate(ordered_hits, start=1):
        ref_start, ref_end = hit_reference_interval(hit)
        rows.append(
            [
                f"<a href=\"#{html.escape(hit_anchor_ids[hit])}\">{html.escape(str(index))}</a>",
                f"<a href=\"#{html.escape(hit_anchor_ids[hit])}\"><code>{html.escape(hit.contig_id)}</code></a>",
                f"<a href=\"#{html.escape(hit_anchor_ids[hit])}-text\">{html.escape(str(hit.contig_length))}</a>",
                f"<a href=\"#{html.escape(hit_anchor_ids[hit])}-text\">{html.escape(str(hit.query_start))}</a>",
                f"<a href=\"#{html.escape(hit_anchor_ids[hit])}-text\">{html.escape(str(hit.query_end))}</a>",
                f"<a href=\"#{html.escape(hit_anchor_ids[hit])}-text\">{html.escape(str(ref_start))}</a>",
                f"<a href=\"#{html.escape(hit_anchor_ids[hit])}-text\">{html.escape(str(ref_end))}</a>",
                f"<a href=\"#{html.escape(hit_anchor_ids[hit])}-text\">{html.escape(describe_hit_strand(hit))}</a>",
                f"<a href=\"#{html.escape(hit_anchor_ids[hit])}-text\">{html.escape(format_report_float(hit.percent_identity, 2))}</a>",
                f"<a href=\"#{html.escape(hit_anchor_ids[hit])}-text\">{html.escape(f'{hit.evalue:.6g}')}</a>",
                f"<a href=\"#{html.escape(hit_anchor_ids[hit])}-text\">{html.escape(format_report_float(hit.query_coverage, 2))}</a>",
                f"<a href=\"#{html.escape(hit_anchor_ids[hit])}\">{html.escape(hit.description)}</a>",
            ]
        )

    body_parts.append("<div class=\"section\" id=\"hit-table\"><h2>Contig Hits</h2></div>")
    body_parts.append(render_reference_detail_hits_table(headers, rows, [detail_hit_row_ids[hit] for hit in ordered_hits]))
    alignment_block_index = render_alignment_block_index(
        ordered_hits,
        contig_anchor_ids,
        hit_anchor_ids,
        detail_hit_row_ids,
        index_row_ids,
    )
    if alignment_block_index:
        body_parts.append(
            "<div class=\"section\" id=\"alignment-index\"><h2>Alignment Index</h2>"
            "<p class=\"muted\">Quick links into the contig alignments below.</p></div>"
        )
        body_parts.append(alignment_block_index)
    alignment_blocks = render_alignment_blocks(
        ordered_hits,
        hit_anchor_ids,
        contig_anchor_ids,
        detail_hit_row_ids,
        index_row_ids,
    )
    if alignment_blocks:
        body_parts.append(alignment_blocks)
    return build_html_document(title, "\n".join(body_parts))


def write_reference_detail_pages(
    output_dir: Path,
    analysis: str,
    summary_rows: list[dict[str, object]],
    hits: list[BlastHit],
) -> str:
    detail_dir = output_dir / f"{analysis}_references"
    detail_dir.mkdir(parents=True, exist_ok=True)
    summary_by_group = {str(row["group_id"]): row for row in summary_rows if row["passes_filters"]}
    grouped_hits: dict[str, list[BlastHit]] = defaultdict(list)
    for hit in hits:
        grouped_hits[hit_group_id(hit)].append(hit)

    for group_id, summary_row in summary_by_group.items():
        group_hits = grouped_hits.get(group_id, [])
        reference_id = str(summary_row["reference_id"])
        page = build_reference_detail_html(
            analysis,
            group_id,
            reference_id,
            group_hits,
            summary_row,
            summary_page_name=f"{analysis}.html",
        )
        write_text_report(detail_dir / f"{group_id}.html", page)

    return detail_dir.name


def build_undetermined_report_html(
    undetermined_records: list[FastaRecord],
    contig_depths: dict[str, ContigDepth],
    best_raw_hits: dict[str, BlastHit],
    read_length_stats: dict[str, dict[int, int]] | None = None,
    sirna_percent: float | None = None,
    hits_only: bool = False,
    reference_targets: dict[tuple[str, str], dict[str, str]] | None = None,
) -> str:
    def build_summary_count_line(label: str, count: int, href: str | None = None) -> str:
        count_text = html.escape(str(count))
        if href is not None:
            count_text = f"<a href=\"{html.escape(href)}\">{count_text}</a>"
        return f"{html.escape(label)}: {count_text}"

    title = "Undetermined Contigs"
    total_count = len(undetermined_records)
    hit_count = sum(1 for record in undetermined_records if record.seq_id in best_raw_hits)
    no_hit_count = total_count - hit_count
    show_sirna_table = bool(read_length_stats)
    threshold_note = ""
    if show_sirna_table and sirna_percent is not None:
        threshold_note = f" Candidate threshold: >= {format_report_float(sirna_percent * 100, 2)}%."
    body_parts = [
        f"<h1 id=\"detail-top\">{title}</h1>",
        (
            "<p class=\"muted\">This summary lists contigs that were not assigned to the final known or novel virus sets. "
            "Candidate contigs are highlighted in green and marked in the Candidate column."
            f"{html.escape(threshold_note)}</p>"
            if show_sirna_table
            else "<p class=\"muted\">This summary lists contigs that were not assigned to the final known or novel virus sets. siRNA size-distribution columns are not available for this run.</p>"
        ),
    ]
    base_headers = [
        "ID",
        "Length",
        "Covered Bases",
        "Contig Coverage",
    ]
    if show_sirna_table:
        base_headers.extend([str(length) for length in range(18, 34)])
        base_headers.append("21-22 (%)")
        base_headers.append("Candidate")
    base_headers.extend(
        [
        "Depth",
        "Depth (Norm)",
        ]
    )
    hit_headers = [
        *base_headers,
        "Acc#",
        "Genus",
        "Description",
        "E value",
        "BLAST",
        "Match ID",
        "Identity",
        "Query Coverage",
    ]
    sortable_rows_with_hits: list[tuple[tuple[float, int, str], str, list[str], str]] = []
    sortable_rows_without_hits: list[tuple[tuple[float, int, str], str, list[str], str]] = []
    candidate_hit_count = 0
    candidate_no_hit_count = 0
    for record in undetermined_records:
        best_hit = best_raw_hits.get(record.seq_id)
        depth = contig_depths.get(record.seq_id, ContigDepth(record.seq_id, len(record.sequence), 0, 0, 0.0, 0.0))
        contig_length = len(record.sequence)
        contig_coverage = (depth.covered_bases / contig_length) if contig_length else 0.0
        base_row = [
            f"<code>{html.escape(record.seq_id)}</code>",
            html.escape(str(contig_length)),
            html.escape(str(depth.covered_bases)),
            html.escape(f"{format_report_float(contig_coverage * 100, 1)}%"),
        ]
        row_class = ""
        is_candidate = False
        if show_sirna_table:
            length_counts = (read_length_stats or {}).get(record.seq_id, {})
            total_sirna = sum(length_counts.get(length, 0) for length in range(18, 34))
            ratio_value = (
                (length_counts.get(21, 0) + length_counts.get(22, 0)) / total_sirna
                if total_sirna
                else None
            )
            base_row.extend(html.escape(str(length_counts.get(length, 0))) for length in range(18, 34))
            base_row.append(
                html.escape(f"{format_report_float(ratio_value * 100, 2)}%") if ratio_value is not None else "NA"
            )
            if ratio_value is not None and sirna_percent is not None and ratio_value >= sirna_percent:
                is_candidate = True
                row_class = "candidate-row"
            base_row.append(
                "<span class=\"candidate-flag candidate-yes\">Yes</span>"
                if is_candidate
                else "<span class=\"candidate-flag candidate-no\">No</span>"
            )
        base_row = [
            *base_row,
            html.escape(format_report_float(depth.mean_depth, 2)),
            html.escape(format_report_float(depth.normalized_depth, 2)),
        ]
        sort_key = (-depth.normalized_depth, -contig_length, record.seq_id)
        if best_hit is not None:
            reference_target = (reference_targets or {}).get((best_hit.analysis, hit_group_id(best_hit)))
            best_hit_text = f"<code>{html.escape(best_hit.hit_id)}</code>"
            reference_id_text = f"<code>{html.escape(best_hit.reference_id)}</code>"
            analysis_text = html.escape(best_hit.analysis)
            description_text = html.escape(best_hit.description)
            if reference_target is not None:
                detail_href = reference_target.get("detail_href")
                summary_href = reference_target.get("summary_href")
                if detail_href:
                    best_hit_text = f"<a href=\"{html.escape(detail_href)}\">{best_hit_text}</a>"
                    reference_id_text = f"<a href=\"{html.escape(detail_href)}\">{reference_id_text}</a>"
                    description_text = f"<a href=\"{html.escape(detail_href)}\">{description_text}</a>"
                if summary_href:
                    analysis_text = f"<a href=\"{html.escape(summary_href)}\">{analysis_text}</a>"
            if is_candidate:
                candidate_hit_count += 1
            sortable_rows_with_hits.append(
                (
                    sort_key,
                    record.seq_id,
                    [
                        *base_row,
                        reference_id_text,
                        html.escape(best_hit.genus),
                        description_text,
                        html.escape(f"{best_hit.evalue:.6g}"),
                        analysis_text,
                        best_hit_text,
                        html.escape(format_report_float(best_hit.percent_identity, 2)),
                        html.escape(format_report_float(best_hit.query_coverage, 2)),
                    ],
                    row_class,
                )
            )
        else:
            if is_candidate:
                candidate_no_hit_count += 1
            sortable_rows_without_hits.append((sort_key, record.seq_id, base_row, row_class))

    sorted_rows_with_hits = sorted(sortable_rows_with_hits)
    sorted_rows_without_hits = sorted(sortable_rows_without_hits)

    contig_ids_with_hits = [contig_id for _, contig_id, _, _ in sorted_rows_with_hits]
    contig_ids_without_hits = [contig_id for _, contig_id, _, _ in sorted_rows_without_hits]
    row_ids_with_hits = [build_undetermined_row_id("hit", contig_id) for contig_id in contig_ids_with_hits]
    row_ids_without_hits = [build_undetermined_row_id("no-match", contig_id) for contig_id in contig_ids_without_hits]

    rows_with_hits = []
    row_classes_with_hits = []
    for (_, contig_id, row, row_class), row_id in zip(sorted_rows_with_hits, row_ids_with_hits):
        linked_row = row.copy()
        linked_row[0] = build_undetermined_contig_cell(contig_id, row_id)
        rows_with_hits.append(linked_row)
        row_classes_with_hits.append(row_class)

    rows_without_hits = []
    row_classes_without_hits = []
    for (_, contig_id, row, row_class), row_id in zip(sorted_rows_without_hits, row_ids_without_hits):
        linked_row = row.copy()
        linked_row[0] = build_undetermined_contig_cell(contig_id, row_id)
        rows_without_hits.append(linked_row)
        row_classes_without_hits.append(row_class)

    candidate_hit_indices = [index for index, row_class in enumerate(row_classes_with_hits) if row_class == "candidate-row"]
    candidate_no_hit_indices = [
        index for index, row_class in enumerate(row_classes_without_hits) if row_class == "candidate-row"
    ]
    candidate_rows_with_hits = []
    candidate_row_ids_with_hits = []
    candidate_row_classes_with_hits = []
    for index in candidate_hit_indices:
        candidate_row = rows_with_hits[index].copy()
        candidate_row[0] = build_undetermined_contig_cell(contig_ids_with_hits[index], row_ids_with_hits[index])
        candidate_rows_with_hits.append(candidate_row)
        candidate_row_ids_with_hits.append(build_undetermined_row_id("candidate-hit", contig_ids_with_hits[index]))
        candidate_row_classes_with_hits.append(row_classes_with_hits[index])

    candidate_rows_without_hits = []
    candidate_row_ids_without_hits = []
    candidate_row_classes_without_hits = []
    for index in candidate_no_hit_indices:
        candidate_row = rows_without_hits[index].copy()
        candidate_row[0] = build_undetermined_contig_cell(contig_ids_without_hits[index], row_ids_without_hits[index])
        candidate_rows_without_hits.append(candidate_row)
        candidate_row_ids_without_hits.append(build_undetermined_row_id("candidate-no-match", contig_ids_without_hits[index]))
        candidate_row_classes_without_hits.append(row_classes_without_hits[index])
    candidate_total_count = candidate_hit_count + candidate_no_hit_count

    navigation_links: list[str] = []
    if hits_only:
        if show_sirna_table and candidate_rows_with_hits:
            navigation_links.append("<a href=\"#candidate-table\">Candidate Table</a>")
        navigation_links.append("<a href=\"#hit-table\">Match Table</a>")
        navigation_links.append("<a href=\"undetermined.html#no-match-table\">No-Match View</a>")
    else:
        if show_sirna_table and candidate_rows_without_hits:
            navigation_links.append("<a href=\"#candidate-table\">Candidate Table</a>")
        if no_hit_count:
            navigation_links.append("<a href=\"#no-match-table\">No-Match Table</a>")
        if hit_count:
            navigation_links.append("<a href=\"undetermined_blast.html#hit-table\">Match Table</a>")
    if navigation_links:
        navigation_links.append("<a href=\"#detail-top\">Top</a>")
    if navigation_links:
        body_parts.append("<p class=\"section-nav\">" + " | ".join(navigation_links) + "</p>")

    if not hits_only:
        summary_lines = [
            build_summary_count_line("Total undetermined contigs", total_count),
            build_summary_count_line(
                "With virus-database matches",
                hit_count,
                "undetermined_blast.html#hit-table" if hit_count else None,
            ),
            build_summary_count_line("Without virus-database matches", no_hit_count, "#no-match-table" if no_hit_count else None),
        ]
        if show_sirna_table:
            summary_lines.append(
                build_summary_count_line("21-22 nt candidate rows", candidate_total_count, "#candidate-table" if candidate_total_count else None)
            )
        summary_lines.append("This page shows contigs without virus-database matches; use <code>undetermined_blast.html</code> for contigs with matches.")
        body_parts.append("<p>" + "<br>".join(summary_lines) + "</p>")
    else:
        summary_lines = [build_summary_count_line("Undetermined contigs with virus-database matches", hit_count, "#hit-table" if hit_count else None)]
        if show_sirna_table:
            summary_lines.append(
                build_summary_count_line(
                    "21-22 nt candidate rows in this view",
                    candidate_hit_count,
                    "#candidate-table" if candidate_hit_count else None,
                )
            )
        body_parts.append("<p>" + "<br>".join(summary_lines) + "</p>")

    if hits_only:
        if not rows_with_hits:
            body_parts.append("<p>No undetermined contigs matched this page.</p>")
            return build_html_document(title, "\n".join(body_parts))
        if show_sirna_table and candidate_rows_with_hits:
            body_parts.append(
                (
                    "<div class=\"section\" id=\"candidate-table\"><h2>Candidate Contigs</h2>"
                    "<p class=\"muted\">Subset of undetermined contigs in this hit view that pass the configured 21-22 nt enrichment threshold. "
                    "Contig IDs jump to the matching row in the full match table below.</p></div>"
                )
            )
            body_parts.append(
                render_undetermined_table(
                    hit_headers,
                    candidate_rows_with_hits,
                    candidate_row_classes_with_hits,
                    candidate_row_ids_with_hits,
                    show_sirna_table=show_sirna_table,
                    include_hit_columns=True,
                )
            )
        body_parts.append(
            (
                "<div class=\"section\" id=\"hit-table\"><h2>With Virus-Database Matches</h2>"
                "<p class=\"muted\">Undetermined contigs with virus-database matches.</p></div>"
            )
        )
        body_parts.append(
            render_undetermined_table(
                hit_headers,
                rows_with_hits,
                row_classes_with_hits,
                row_ids_with_hits,
                show_sirna_table=show_sirna_table,
                include_hit_columns=True,
            )
        )
        return build_html_document(title, "\n".join(body_parts))

    if show_sirna_table and candidate_rows_without_hits:
        body_parts.append(
            (
                "<div class=\"section\" id=\"candidate-table\"><h2>Candidate Contigs</h2>"
                "<p class=\"muted\">Subset of undetermined contigs in this no-match view that pass the configured 21-22 nt enrichment threshold. "
                "Contig IDs jump to the matching row in the full no-match table below.</p></div>"
            )
        )
        body_parts.append(
            render_undetermined_table(
                base_headers,
                candidate_rows_without_hits,
                candidate_row_classes_without_hits,
                candidate_row_ids_without_hits,
                show_sirna_table=show_sirna_table,
                include_hit_columns=False,
            )
        )

    if rows_without_hits:
        body_parts.append(
            (
                "<div class=\"section\" id=\"no-match-table\"><h2>Without Virus-Database Matches</h2>"
                "<p class=\"muted\">Undetermined contigs without virus-database matches.</p></div>"
            )
        )
        body_parts.append(
            render_undetermined_table(
                base_headers,
                rows_without_hits,
                row_classes_without_hits,
                row_ids_without_hits,
                show_sirna_table=show_sirna_table,
                include_hit_columns=False,
            )
        )

    if not rows_without_hits:
        body_parts.append("<p>No undetermined contigs matched this page.</p>")
        return build_html_document(title, "\n".join(body_parts))

    return build_html_document(title, "\n".join(body_parts))


def write_text_report(path: Path, content: str) -> None:
    path.write_text(content, encoding="utf-8")


def write_marker(path: Path, should_exist: bool) -> None:
    if should_exist:
        path.touch()
    elif path.exists():
        path.unlink()


def run_python_identifier(args, db_path: str, sample_path: Path, contig_path: Path, output_dir: Path, data_type: str) -> IdentifyResult:
    result_dir = output_dir / f"result_{sample_path.name}"
    if result_dir.exists():
        shutil.rmtree(result_dir)
    result_dir.mkdir(parents=True, exist_ok=True)

    contig_records = read_fasta_records(contig_path)
    records_by_id = {record.seq_id: record for record in contig_records}
    write_fasta_records(result_dir / "contig_sequences.fa", contig_records)

    tool_paths = tool_path_map()
    nucleotide_reference = resolve_reference_path(args.reference, db_path)
    protein_reference = resolve_protein_reference_path(args.reference, db_path)
    ensure_blast_db(nucleotide_reference, "nucl", tool_paths, debug=args.debug)
    ensure_blast_db(protein_reference, "prot", tool_paths, debug=args.debug)

    annotations, protein_to_reference = load_database_annotations(db_path)
    contig_depths, contig_depth_tsv, read_length_stats = compute_contig_depths(
        sample_path,
        contig_path,
        contig_records,
        result_dir,
        args,
        tool_paths,
    )
    library_size = count_sequences(sample_path)

    blastn_output = result_file(result_dir, "blastn.outfmt6.tsv")
    run_blastn(contig_path, nucleotide_reference, blastn_output, args, tool_paths)
    blastn_hits = parse_blast_hits("blastn", blastn_output, annotations, protein_to_reference)
    blastn_candidate_hits = filter_hits(blastn_hits, BLASTN_MIN_IDENTITY, BLASTN_MIN_QUERY_COVERAGE)
    blastn_top_hits = select_best_hits(blastn_hits, BLASTN_MIN_IDENTITY, BLASTN_MIN_QUERY_COVERAGE)
    blastn_summary_rows = summarize_hits(
        blastn_candidate_hits,
        contig_depths=contig_depths,
        library_size=library_size,
        coverage_cutoff=args.coverage_cutoff,
        depth_cutoff=args.depth_cutoff,
        norm_depth_cutoff=args.norm_depth_cutoff,
    )
    blastn_kept_groups, _ = prune_redundant_group_ids(
        blastn_summary_rows,
        build_group_support(blastn_candidate_hits),
        diff_ratio=args.diff_ratio,
        diff_contig_cover=args.diff_contig_cover,
        diff_contig_length=args.diff_contig_length,
    )
    blastn_summary_rows = apply_summary_filter_flags(blastn_summary_rows, blastn_kept_groups)
    blastn_allowed_groups = {row["group_id"] for row in blastn_summary_rows if row["passes_filters"]}
    blastn_final_hits = select_best_hits(
        blastn_hits,
        BLASTN_MIN_IDENTITY,
        BLASTN_MIN_QUERY_COVERAGE,
        allowed_group_ids=blastn_allowed_groups,
    )
    blastn_reference_records = filter_reference_records(read_fasta_records(nucleotide_reference), blastn_allowed_groups)
    blastn_reference_path = result_dir / "blastn.references.fa"
    write_fasta_records(blastn_reference_path, blastn_reference_records)

    known_ids = {hit.contig_id for hit in blastn_final_hits}
    known_records = filter_records(contig_records, known_ids)
    remaining_records = [record for record in contig_records if record.seq_id not in known_ids]

    blastx_hits: list[BlastHit] = []
    blastx_candidate_hits: list[BlastHit] = []
    blastx_top_hits: list[BlastHit] = []
    blastx_summary_rows: list[dict[str, object]] = []
    blastx_final_hits: list[BlastHit] = []
    blastx_allowed_groups: set[str] = set()
    blastx_reference_path = result_dir / "blastx.references.fa"
    write_fasta_records(blastx_reference_path, [])
    remaining_path = result_file(result_dir, "undetermined.input.fa")
    if remaining_records:
        write_fasta_records(remaining_path, remaining_records)
        blastx_output = result_file(result_dir, "blastx.outfmt6.tsv")
        effective_exp_valuex = args.exp_valuex if args.exp_valuex is not None else (1e-5 if data_type == "mRNA" else 1e-2)
        run_blastx(remaining_path, protein_reference, blastx_output, args, effective_exp_valuex, tool_paths)
        blastx_hits = parse_blast_hits("blastx", blastx_output, annotations, protein_to_reference)
        blastx_candidate_hits = filter_hits(blastx_hits, args.percent_identity, BLASTX_MIN_QUERY_COVERAGE)
        blastx_top_hits = select_best_hits(blastx_hits, args.percent_identity, BLASTX_MIN_QUERY_COVERAGE)
        blastx_summary_rows = summarize_hits(
            blastx_candidate_hits,
            contig_depths=contig_depths,
            library_size=library_size,
            coverage_cutoff=args.coverage_cutoff,
            depth_cutoff=args.depth_cutoff,
            norm_depth_cutoff=args.norm_depth_cutoff,
        )
        blastx_kept_groups, _ = prune_redundant_group_ids(
            blastx_summary_rows,
            build_group_support(blastx_candidate_hits),
            diff_ratio=args.diff_ratio,
            diff_contig_cover=args.diff_contig_cover,
            diff_contig_length=args.diff_contig_length,
        )
        blastx_summary_rows = apply_summary_filter_flags(blastx_summary_rows, blastx_kept_groups)
        blastx_allowed_groups = {row["group_id"] for row in blastx_summary_rows if row["passes_filters"]}
        blastx_final_hits = select_best_hits(
            blastx_hits,
            args.percent_identity,
            BLASTX_MIN_QUERY_COVERAGE,
            allowed_group_ids=blastx_allowed_groups,
        )
        blastx_reference_records = filter_reference_records(read_fasta_records(protein_reference), blastx_allowed_groups)
        write_fasta_records(blastx_reference_path, blastx_reference_records)
    else:
        remaining_path.write_text("", encoding="utf-8")

    novel_ids = {hit.contig_id for hit in blastx_final_hits}
    novel_records = [record for record in remaining_records if record.seq_id in novel_ids]
    undetermined_records = [record for record in remaining_records if record.seq_id not in novel_ids]

    write_fasta_records(result_dir / "contig_sequences.blastn.fa", known_records)
    write_fasta_records(result_dir / "contig_sequences.blastx.fa", novel_records)
    write_fasta_records(result_dir / "contig_sequences.undetermined.fa", undetermined_records)

    blastn_report_hits = [
        hit for hit in blastn_candidate_hits if hit_group_id(hit) in blastn_allowed_groups and hit.contig_id in known_ids
    ]
    blastx_report_hits = [hit for hit in blastx_candidate_hits if hit_group_id(hit) in blastx_allowed_groups]
    write_legacy_alignment_table(result_file(result_dir, "blastn.xls"), blastn_report_hits, records_by_id)
    write_legacy_alignment_table(result_file(result_dir, "blastx.xls"), blastx_report_hits, records_by_id)
    write_legacy_alignment_sam(result_file(result_dir, "blastn.sam"), blastn_report_hits)
    write_legacy_alignment_sam(result_file(result_dir, "blastx.sam"), blastx_report_hits)

    blastn_raw_tsv = result_file(result_dir, "blastn.raw.tsv")
    blastn_top_tsv = result_file(result_dir, "blastn.top.tsv")
    blastx_raw_tsv = result_file(result_dir, "blastx.raw.tsv")
    blastx_top_tsv = result_file(result_dir, "blastx.top.tsv")
    write_hit_table(blastn_raw_tsv, blastn_hits)
    write_hit_table(blastn_top_tsv, blastn_top_hits)
    write_hit_table(blastx_raw_tsv, blastx_hits)
    write_hit_table(blastx_top_tsv, blastx_top_hits)
    reference_targets = build_reference_target_map(
        ("blastn", blastn_summary_rows),
        ("blastx", blastx_summary_rows),
    )

    blastn_html = result_dir / "blastn.html"
    blastx_html = result_dir / "blastx.html"
    undetermined_html = result_dir / "undetermined.html"
    undetermined_blast_html = result_dir / "undetermined_blast.html"
    write_text_report(
        blastn_html,
        build_reference_report_html(
            "blastn",
            blastn_summary_rows,
            blastn_final_hits,
            blastn_reference_path.name,
            None,
        ),
    )
    write_text_report(
        blastx_html,
        build_reference_report_html(
            "blastx",
            blastx_summary_rows,
            blastx_final_hits,
            blastx_reference_path.name,
            None,
        ),
    )
    best_undetermined_hits = select_best_raw_hits(
        [hit for hit in blastn_hits + blastx_hits if hit.contig_id in {record.seq_id for record in undetermined_records}]
    )
    write_text_report(
        undetermined_html,
        build_undetermined_report_html(
            undetermined_records,
            contig_depths,
            best_undetermined_hits,
            read_length_stats=read_length_stats if data_type == "sRNA" else None,
            sirna_percent=args.siRNA_percent if data_type == "sRNA" else None,
            hits_only=False,
            reference_targets=reference_targets,
        ),
    )
    write_text_report(
        undetermined_blast_html,
        build_undetermined_report_html(
            undetermined_records,
            contig_depths,
            best_undetermined_hits,
            read_length_stats=read_length_stats if data_type == "sRNA" else None,
            sirna_percent=args.siRNA_percent if data_type == "sRNA" else None,
            hits_only=True,
            reference_targets=reference_targets,
        ),
    )

    summary_rows = blastn_summary_rows + blastx_summary_rows
    summary_tsv = ""
    summary_json = ""

    write_marker(result_dir / "no_virus_detected_by_blastn", len(blastn_final_hits) == 0)
    write_marker(result_dir / "no_virus_detected_by_blastx", len(blastx_final_hits) == 0)

    cleanup_result_dir(
        result_dir,
        keep_names={
            "contig_sequences.fa",
            "contig_sequences.blastn.fa",
            "contig_sequences.blastx.fa",
            "contig_sequences.undetermined.fa",
            "blastn.references.fa",
            "blastx.references.fa",
            "blastn.html",
            "blastx.html",
            "undetermined.html",
            "undetermined_blast.html",
            "blastn.sam",
            "blastx.sam",
            "blastn.xls",
            "blastx.xls",
        },
    )

    return IdentifyResult(
        result_dir=str(result_dir),
        blastn_raw_tsv=str(blastn_raw_tsv),
        blastn_top_tsv=str(blastn_top_tsv),
        blastx_raw_tsv=str(blastx_raw_tsv),
        blastx_top_tsv=str(blastx_top_tsv),
        summary_tsv=str(summary_tsv),
        summary_json=str(summary_json),
        known_contig_count=len(known_records),
        novel_contig_count=len(novel_records),
        undetermined_contig_count=len(undetermined_records),
    )
