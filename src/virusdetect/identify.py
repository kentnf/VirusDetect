from __future__ import annotations

import csv
import gzip
import html
import json
import subprocess
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path

from virusdetect.db import resolve_data_path, resolve_reference_path, resolve_seq_info_path
from virusdetect.fasta import FastaRecord, read_fasta_records, write_fasta_records
from virusdetect.runtime import tool_path_map
from virusdetect.sample import count_sequences
from virusdetect.sam import expand_xa_hits

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
)
BLASTN_MIN_IDENTITY = 60.0
BLASTN_MIN_QUERY_COVERAGE = 50.0
BLASTX_MIN_QUERY_COVERAGE = 60.0


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


def compute_contig_depths(
    sample_path: Path,
    contig_path: Path,
    contig_records: list[FastaRecord],
    result_dir: Path,
    args,
    tool_paths: dict[str, str],
) -> tuple[dict[str, ContigDepth], Path]:
    depth_path = result_dir / f"{sample_path.name}.contig_depth.tsv"
    if not contig_records:
        write_contig_depth_tsv(depth_path, {})
        return {}, depth_path

    bwa = tool_paths.get("bwa")
    samtools = tool_paths.get("samtools")
    if not bwa or not samtools:
        raise SystemExit("bwa and samtools are required for contig-depth estimation")

    reference_path = str(contig_path)
    ensure_bwa_index(reference_path, {"bwa": bwa}, debug=args.debug)
    ensure_fasta_index(reference_path, {"samtools": samtools}, debug=args.debug)

    sai_path = result_dir / f"{sample_path.name}.contig_depth.bwa.sai"
    sam_path = result_dir / f"{sample_path.name}.contig_depth.sam"
    bam_path = result_dir / f"{sample_path.name}.contig_depth.bam"
    sorted_prefix = result_dir / f"{sample_path.name}.contig_depth.sorted"
    sorted_bam_path = result_dir / f"{sample_path.name}.contig_depth.sorted.bam"
    raw_depth_path = result_dir / f"{sample_path.name}.contig_depth.pileup"

    align_command = build_bwa_align_command(reference_path, str(sample_path), args.threads)
    align_command[0] = bwa
    run_command_to_file(align_command, sai_path, debug=args.debug)

    samse_command = build_bwa_samse_command(reference_path, str(sai_path), str(sample_path))
    samse_command[0] = bwa
    run_command_to_file(samse_command, sam_path, debug=args.debug)
    expand_xa_hits(sam_path)

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
    return contig_depths, depth_path


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
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, field_names, delimiter="\t")
        writer.writeheader()
        for hit in hits:
            writer.writerow(asdict(hit))


def write_legacy_style_table(path: Path, hits: list[BlastHit], records_by_id: dict[str, FastaRecord]) -> None:
    field_names = [
        "Contig_ID",
        "Contig_Seq",
        "Contig_Len",
        "Hit_ID",
        "Reference_ID",
        "Hit_Len",
        "Genus",
        "Description",
        "Contig_start",
        "Contig_end",
        "Hit_start",
        "Hit_end",
        "Hsp_identity",
        "E_value",
        "Query_coverage",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, field_names, delimiter="\t")
        writer.writeheader()
        for hit in hits:
            record = records_by_id.get(hit.contig_id)
            writer.writerow(
                {
                    "Contig_ID": hit.contig_id,
                    "Contig_Seq": record.sequence if record else "",
                    "Contig_Len": hit.contig_length,
                    "Hit_ID": hit.hit_id,
                    "Reference_ID": hit.reference_id,
                    "Hit_Len": hit.hit_length,
                    "Genus": hit.genus,
                    "Description": hit.description,
                    "Contig_start": hit.query_start,
                    "Contig_end": hit.query_end,
                    "Hit_start": hit.hit_start,
                    "Hit_end": hit.hit_end,
                    "Hsp_identity": f"{hit.percent_identity:.2f}",
                    "E_value": f"{hit.evalue:.6g}",
                    "Query_coverage": f"{hit.query_coverage:.2f}",
                }
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


def render_coverage_blocks(reference_length: int, hits: list[BlastHit]) -> str:
    blocks: list[str] = []
    for hit in sorted(hits, key=lambda item: (hit_reference_interval(item)[0], hit_reference_interval(item)[1], item.contig_id)):
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


def render_reference_coverage_map(reference_length: int, hits: list[BlastHit]) -> str:
    if reference_length <= 0 or not hits:
        return "<p class=\"muted\">Reference coverage map is unavailable for this reference.</p>"

    return (
        "<div class=\"section\">"
        "<h2>Reference Coverage Map</h2>"
        f"<div class=\"coverage-axis\"><span>1</span><span>{reference_length}</span></div>"
        f"<div class=\"coverage-track\">{render_coverage_blocks(reference_length, hits)}</div>"
        "<p class=\"muted\">Each block shows one contig aligned on the reference coordinate axis. Blue indicates forward strand; orange indicates reverse strand.</p>"
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
        "    .coverage-axis { display: flex; justify-content: space-between; font-size: 12px; color: #52606d; margin: 6px 0; }\n"
        "    .coverage-track { position: relative; height: 20px; border-radius: 999px; background: #d9e2ec; overflow: hidden; }\n"
        "    .coverage-track-compact { min-width: 180px; height: 14px; }\n"
        "    .coverage-block { position: absolute; top: 0; height: 100%; border-radius: 999px; opacity: 0.9; }\n"
        "    .coverage-block.forward { background: #2f6fed; }\n"
        "    .coverage-block.reverse { background: #d97706; }\n"
        "  </style>\n"
        "</head>\n"
        "<body>\n"
        f"{body}\n"
        "</body>\n"
        "</html>\n"
    )


def render_html_table(headers: list[str], rows: list[list[str]]) -> str:
    header_html = "".join(f"<th>{html.escape(header)}</th>" for header in headers)
    body_rows = []
    for row in rows:
        body_rows.append("<tr>" + "".join(f"<td>{cell}</td>" for cell in row) + "</tr>")
    return "<table>\n<thead><tr>" + header_html + "</tr></thead>\n<tbody>\n" + "\n".join(body_rows) + "\n</tbody>\n</table>"


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
        f"<h1>{html.escape(title)}</h1>",
        "<p class=\"muted\">Minimal Python-native summary report. Detailed per-reference graphics from the Perl workflow are not yet ported.</p>",
        f"<p>Reference FASTA: <code>{html.escape(reference_fasta_name)}</code></p>",
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
    for row in kept_rows:
        group_id = str(row["group_id"])
        stats = identity_stats.get(group_id, {})
        group_hits = grouped_hits.get(group_id, [])
        rows.append(
            [
                (
                    f"<a href=\"{html.escape(detail_dir_name)}/{html.escape(group_id)}.html\"><code>{html.escape(group_id)}</code></a>"
                    if detail_dir_name
                    else f"<code>{html.escape(group_id)}</code>"
                ),
                f"<code>{html.escape(str(row['reference_id']))}</code>",
                html.escape(str(row["hit_length"])),
                html.escape(f"{row['covered_bases']} ({float(row['reference_coverage']) * 100:.1f}%)"),
                render_reference_coverage_strip(int(row["hit_length"]), group_hits),
                html.escape(str(row["contig_count"])),
                html.escape(format_report_float(float(row["reference_depth"]), 1)),
                html.escape(format_report_float(float(row["normalized_depth"]), 1)),
                html.escape(format_report_float(float(stats.get("mean", row["mean_percent_identity"])), 2)),
                html.escape(format_report_float(float(stats.get("max", row["mean_percent_identity"])), 2)),
                html.escape(format_report_float(float(stats.get("min", row["mean_percent_identity"])), 2)),
                html.escape(str(row["genus"])),
                html.escape(str(row["description"])),
            ]
        )

    body_parts.append(render_html_table(headers, rows))
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
    body_parts = [
        f"<h1>{html.escape(title)}</h1>",
        f"<p class=\"muted\">Reference detail page for <code>{html.escape(group_id)}</code>.</p>",
    ]

    if summary_row is not None:
        reference_length = int(summary_row.get("hit_length", 0))
        body_parts.append(
            "<p>"
            f"Reference ID: <code>{html.escape(reference_id)}</code><br>"
            f"Reference Length: {html.escape(str(reference_length))}<br>"
            f"Coverage: {html.escape(str(summary_row['covered_bases']))} bases "
            f"({float(summary_row['reference_coverage']) * 100:.1f}%)<br>"
            f"Contigs: {html.escape(str(summary_row['contig_count']))}<br>"
            f"Depth: {html.escape(format_report_float(float(summary_row['reference_depth']), 1))}<br>"
            f"Depth (Norm): {html.escape(format_report_float(float(summary_row['normalized_depth']), 1))}"
            "</p>"
        )
    if summary_page_name:
        body_parts.append(
            f"<p><a href=\"../{html.escape(summary_page_name)}\">Back to {html.escape(analysis.upper())} summary</a></p>"
        )

    if not hits:
        body_parts.append("<p>No contig hits were recorded for this reference.</p>")
        return build_html_document(title, "\n".join(body_parts))

    if summary_row is not None:
        body_parts.append(render_reference_coverage_map(int(summary_row.get("hit_length", 0)), ordered_hits))

    headers = [
        "Order",
        "Contig ID",
        "Contig Length",
        "Query Start",
        "Query End",
        "Ref Range",
        "Ref Span",
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
                html.escape(str(index)),
                f"<code>{html.escape(hit.contig_id)}</code>",
                html.escape(str(hit.contig_length)),
                html.escape(str(hit.query_start)),
                html.escape(str(hit.query_end)),
                html.escape(f"{ref_start}-{ref_end}"),
                html.escape(str(ref_end - ref_start + 1)),
                html.escape(describe_hit_strand(hit)),
                html.escape(format_report_float(hit.percent_identity, 2)),
                html.escape(f"{hit.evalue:.6g}"),
                html.escape(format_report_float(hit.query_coverage, 2)),
                html.escape(hit.description),
            ]
        )

    body_parts.append(render_html_table(headers, rows))
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
    hits_only: bool = False,
) -> str:
    title = "Undetermined Contigs"
    total_count = len(undetermined_records)
    hit_count = sum(1 for record in undetermined_records if record.seq_id in best_raw_hits)
    no_hit_count = total_count - hit_count
    body_parts = [
        f"<h1>{title}</h1>",
        "<p class=\"muted\">This Python-native summary lists contigs that were not assigned to the final known or novel virus sets. siRNA size-distribution reporting is not yet ported.</p>",
        (
            f"<p>Total undetermined contigs: {total_count}<br>"
            f"With virus-database hits: {hit_count}<br>"
            f"Without virus-database hits: {no_hit_count}</p>"
            if not hits_only
            else f"<p>Undetermined contigs with virus-database hits: {hit_count}</p>"
        ),
    ]

    headers = [
        "Contig ID",
        "Length",
        "Covered Bases",
        "Contig Coverage",
        "Mean Depth",
        "Depth (Norm)",
        "Best Hit",
        "Reference ID",
        "Analysis",
        "E-value",
        "Identity",
        "Query Coverage",
        "Description",
    ]
    rows_with_hits: list[list[str]] = []
    rows_without_hits: list[list[str]] = []
    sortable_rows_with_hits: list[tuple[tuple[float, int, str], list[str]]] = []
    sortable_rows_without_hits: list[tuple[tuple[float, int, str], list[str]]] = []
    for record in undetermined_records:
        best_hit = best_raw_hits.get(record.seq_id)
        depth = contig_depths.get(record.seq_id, ContigDepth(record.seq_id, len(record.sequence), 0, 0, 0.0, 0.0))
        contig_length = len(record.sequence)
        contig_coverage = (depth.covered_bases / contig_length) if contig_length else 0.0
        row = [
            f"<code>{html.escape(record.seq_id)}</code>",
            html.escape(str(contig_length)),
            html.escape(str(depth.covered_bases)),
            html.escape(f"{format_report_float(contig_coverage * 100, 1)}%"),
            html.escape(format_report_float(depth.mean_depth, 2)),
            html.escape(format_report_float(depth.normalized_depth, 2)),
            f"<code>{html.escape(best_hit.hit_id)}</code>" if best_hit else "-",
            f"<code>{html.escape(best_hit.reference_id)}</code>" if best_hit else "-",
            html.escape(best_hit.analysis) if best_hit else "-",
            html.escape(f"{best_hit.evalue:.6g}") if best_hit else "-",
            html.escape(format_report_float(best_hit.percent_identity, 2)) if best_hit else "-",
            html.escape(format_report_float(best_hit.query_coverage, 2)) if best_hit else "-",
            html.escape(best_hit.description) if best_hit else "-",
        ]
        sort_key = (-depth.normalized_depth, -contig_length, record.seq_id)
        if best_hit is not None:
            sortable_rows_with_hits.append((sort_key, row))
        else:
            sortable_rows_without_hits.append((sort_key, row))

    rows_with_hits = [row for _, row in sorted(sortable_rows_with_hits)]
    rows_without_hits = [row for _, row in sorted(sortable_rows_without_hits)]

    if hits_only:
        if not rows_with_hits:
            body_parts.append("<p>No undetermined contigs matched this report view.</p>")
            return build_html_document(title, "\n".join(body_parts))
        body_parts.append("<div class=\"section\"><h2>With Virus-Database Hits</h2></div>")
        body_parts.append(render_html_table(headers, rows_with_hits))
        return build_html_document(title, "\n".join(body_parts))

    if rows_with_hits:
        body_parts.append(
            "<div class=\"section\"><h2>With Virus-Database Hits</h2><p class=\"muted\">These contigs matched viral references but were not retained in the final known or novel sets.</p></div>"
        )
        body_parts.append(render_html_table(headers, rows_with_hits))
    if rows_without_hits:
        body_parts.append(
            "<div class=\"section\"><h2>Without Virus-Database Hits</h2><p class=\"muted\">These contigs remained undetermined with no retained best hit in the virus databases.</p></div>"
        )
        body_parts.append(render_html_table(headers, rows_without_hits))

    if not rows_with_hits and not rows_without_hits:
        body_parts.append("<p>No undetermined contigs matched this report view.</p>")
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
    contig_depths, contig_depth_tsv = compute_contig_depths(
        sample_path,
        contig_path,
        contig_records,
        result_dir,
        args,
        tool_paths,
    )
    library_size = count_sequences(sample_path)

    blastn_output = result_dir / f"{sample_path.name}.blastn.outfmt6.tsv"
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
    blastn_reference_path = result_dir / "blastn.reference.fa"
    write_fasta_records(blastn_reference_path, blastn_reference_records)

    known_ids = {hit.contig_id for hit in blastn_final_hits}
    known_records = filter_records(contig_records, known_ids)
    remaining_records = [record for record in contig_records if record.seq_id not in known_ids]

    blastx_hits: list[BlastHit] = []
    blastx_top_hits: list[BlastHit] = []
    blastx_summary_rows: list[dict[str, object]] = []
    blastx_final_hits: list[BlastHit] = []
    blastx_reference_path = result_dir / "blastx.reference.fa"
    write_fasta_records(blastx_reference_path, [])
    remaining_path = result_dir / f"{sample_path.name}.undetermined.input.fa"
    if remaining_records:
        write_fasta_records(remaining_path, remaining_records)
        blastx_output = result_dir / f"{sample_path.name}.blastx.outfmt6.tsv"
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

    blastn_raw_tsv = result_dir / f"{sample_path.name}.blastn.raw.tsv"
    blastn_top_tsv = result_dir / f"{sample_path.name}.blastn.top.tsv"
    blastx_raw_tsv = result_dir / f"{sample_path.name}.blastx.raw.tsv"
    blastx_top_tsv = result_dir / f"{sample_path.name}.blastx.top.tsv"
    write_hit_table(blastn_raw_tsv, blastn_hits)
    write_hit_table(blastn_top_tsv, blastn_top_hits)
    write_hit_table(blastx_raw_tsv, blastx_hits)
    write_hit_table(blastx_top_tsv, blastx_top_hits)
    write_legacy_style_table(result_dir / f"{sample_path.name}.blastn.xls", blastn_final_hits, records_by_id)
    write_legacy_style_table(result_dir / f"{sample_path.name}.blastx.xls", blastx_final_hits, records_by_id)

    blastn_reference_dir = write_reference_detail_pages(result_dir, "blastn", blastn_summary_rows, blastn_final_hits)
    blastx_reference_dir = write_reference_detail_pages(result_dir, "blastx", blastx_summary_rows, blastx_final_hits)

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
            blastn_reference_dir,
        ),
    )
    write_text_report(
        blastx_html,
        build_reference_report_html(
            "blastx",
            blastx_summary_rows,
            blastx_final_hits,
            blastx_reference_path.name,
            blastx_reference_dir,
        ),
    )
    best_undetermined_hits = select_best_raw_hits(
        [hit for hit in blastn_hits + blastx_hits if hit.contig_id in {record.seq_id for record in undetermined_records}]
    )
    write_text_report(
        undetermined_html,
        build_undetermined_report_html(undetermined_records, contig_depths, best_undetermined_hits, hits_only=False),
    )
    write_text_report(
        undetermined_blast_html,
        build_undetermined_report_html(undetermined_records, contig_depths, best_undetermined_hits, hits_only=True),
    )

    summary_rows = blastn_summary_rows + blastx_summary_rows
    summary_tsv = result_dir / f"{sample_path.name}.summary.tsv"
    summary_json = result_dir / f"{sample_path.name}.summary.json"
    write_summary_tsv(summary_tsv, summary_rows)
    with summary_json.open("w", encoding="utf-8") as handle:
        json.dump(
            {
                "sample": sample_path.name,
                "contig_source": str(contig_path),
                "contig_depth_tsv": str(contig_depth_tsv),
                "blastn_reference_fasta": str(blastn_reference_path),
                "blastx_reference_fasta": str(blastx_reference_path),
                "blastn_reference_dir": str(result_dir / blastn_reference_dir),
                "blastx_reference_dir": str(result_dir / blastx_reference_dir),
                "blastn_html": str(blastn_html),
                "blastx_html": str(blastx_html),
                "undetermined_html": str(undetermined_html),
                "undetermined_blast_html": str(undetermined_blast_html),
                "known_contig_count": len(known_records),
                "novel_contig_count": len(novel_records),
                "undetermined_contig_count": len(undetermined_records),
                "summary": summary_rows,
            },
            handle,
            indent=2,
            sort_keys=True,
        )

    write_marker(result_dir / "no_virus_detected_by_blastn", len(blastn_final_hits) == 0)
    write_marker(result_dir / "no_virus_detected_by_blastx", len(blastx_final_hits) == 0)

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
