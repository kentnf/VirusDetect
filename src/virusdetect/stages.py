from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from virusdetect.db import resolve_data_path, resolve_reference_path, resolve_seq_info_path
from virusdetect.fasta import FastaRecord, read_fasta_records, reverse_complement, write_fasta_records
from virusdetect.pileup import count_fasta_sequences, filter_pileup_by_coverage, pileup_to_contig
from virusdetect.runtime import tool_path_map
from virusdetect.sample import detect_file_type, remove_duplicate_reads
from virusdetect.sam import expand_xa_hits, filter_sam_best_hits, generate_unmapped_reads_from_sam


def run_command(command: list[str], debug: bool = False) -> None:
    if debug:
        print("CMD:", " ".join(command))
    completed = subprocess.run(command, check=False)
    if completed.returncode != 0:
        raise SystemExit(f"Command failed with exit code {completed.returncode}: {' '.join(command)}")


def run_logged_command(command: list[str], log_path: Path, debug: bool = False) -> int:
    if debug:
        print("CMD:", " ".join(command), ">>", log_path)
    with log_path.open("a", encoding="utf-8") as log_handle:
        completed = subprocess.run(command, stdout=log_handle, stderr=subprocess.STDOUT, check=False)
    return completed.returncode


def build_bwa_align_command(tool_paths: dict[str, str], reference_path: str, sample_path: str, sai_path: str, threads: int) -> list[str]:
    return [
        tool_paths["bwa"],
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


def build_bwa_samse_command(tool_paths: dict[str, str], reference_path: str, sai_path: str, sample_path: str, multi_hits: int) -> list[str]:
    return [
        tool_paths["bwa"],
        "samse",
        "-n",
        str(multi_hits),
        reference_path,
        str(sai_path),
        sample_path,
    ]


def run_virus_alignment(sample_path: str, db_path: str, args, temp_dir: Path) -> tuple[Path, int]:
    tool_paths = tool_path_map()
    if "bwa" not in tool_paths:
        raise SystemExit("bwa is required for the Python virus-alignment stage")

    reference_path = resolve_reference_path(args.reference, db_path)
    _seq_info_path = resolve_seq_info_path(db_path)
    temp_dir.mkdir(parents=True, exist_ok=True)
    sample_path_obj = Path(sample_path)
    sam_path = temp_dir / f"{sample_path_obj.name}.virus.sam"
    sai_path = temp_dir / "bwa.sai"

    if not Path(f"{reference_path}.amb").exists():
        run_command([tool_paths["bwa"], "index", "-p", reference_path, "-a", "bwtsw", reference_path], debug=args.debug)

    aln_command = build_bwa_align_command(tool_paths, reference_path, sample_path, str(sai_path), args.threads)
    with sai_path.open("w", encoding="utf-8") as sai_handle:
        if args.debug:
            print("CMD:", " ".join(aln_command), ">", sai_path)
        completed = subprocess.run(aln_command, stdout=sai_handle, check=False)
        if completed.returncode != 0:
            raise SystemExit(f"Command failed with exit code {completed.returncode}: {' '.join(aln_command)}")

    samse_command = build_bwa_samse_command(tool_paths, reference_path, sai_path, sample_path, 10000)
    with sam_path.open("w", encoding="utf-8") as sam_handle:
        if args.debug:
            print("CMD:", " ".join(samse_command), ">", sam_path)
        completed = subprocess.run(samse_command, stdout=sam_handle, check=False)
        if completed.returncode != 0:
            raise SystemExit(f"Command failed with exit code {completed.returncode}: {' '.join(samse_command)}")

    expand_xa_hits(sam_path)
    mapped_num = filter_sam_best_hits(sam_path)
    return sam_path, mapped_num


def ensure_hisat2_index(host_reference_path: str, threads: int, debug: bool = False) -> None:
    host_path = Path(host_reference_path)
    if host_path.with_suffix(".1.ht2").exists() or Path(f"{host_reference_path}.1.ht2").exists():
        return

    tool_paths = tool_path_map()
    hisat2_build = tool_paths.get("hisat2-build") or "hisat2-build"
    run_command([hisat2_build, "-p", str(threads), host_reference_path, host_reference_path], debug=debug)


def run_host_subtraction(sample_path: str, db_path: str, args, temp_dir: Path, file_type: str, data_type: str) -> tuple[Path, int]:
    if not args.host_reference:
        return Path(sample_path), 0

    tool_paths = tool_path_map()
    host_reference_path = resolve_data_path(args.host_reference, db_path)
    sample_name = Path(sample_path).name
    host_sam_path = temp_dir / f"{sample_name}.host.sam"
    unmapped_path = temp_dir / f"{sample_name}.unmapped"

    if data_type == "mRNA":
        ensure_hisat2_index(host_reference_path, args.threads, debug=args.debug)
        hisat_file_type = "-q" if file_type == "fastq" else "-f"
        report_path = temp_dir / f"{sample_name}.hisat.report.txt"
        hisat_command = [
            tool_paths["hisat2"],
            "--time",
            "-p",
            str(args.threads),
            "--un",
            str(unmapped_path),
            "--no-unal",
            hisat_file_type,
            "-k",
            "1",
            "--mp",
            "1,1",
            "--rdg",
            "0,1",
            "--rfg",
            "0,1",
            "--np",
            "1",
            "--score-min",
            f"C,-{args.hisat_dist},0",
            "--ignore-quals",
            "-x",
            host_reference_path,
            "-U",
            sample_path,
            "-S",
            str(host_sam_path),
        ]
        with report_path.open("w", encoding="utf-8") as report_handle:
            if args.debug:
                print("CMD:", " ".join(hisat_command), ">", report_path)
            completed = subprocess.run(hisat_command, stdout=report_handle, stderr=subprocess.STDOUT, check=False)
            if completed.returncode != 0:
                raise SystemExit(f"Command failed with exit code {completed.returncode}: {' '.join(hisat_command)}")

        report = report_path.read_text(encoding="utf-8")
        total_num = 0
        unmap_num = 0
        import re
        total_match = re.search(r"(\d+) reads; of these:", report)
        unmap_match = re.search(r"\s+(\d+) .*aligned 0 times", report)
        if total_match:
            total_num = int(total_match.group(1))
        if unmap_match:
            unmap_num = int(unmap_match.group(1))
        mapped_num = max(total_num - unmap_num, 0)
        return unmapped_path, mapped_num

    if "bwa" not in tool_paths:
        raise SystemExit("bwa is required for sRNA host subtraction")

    if not Path(f"{host_reference_path}.amb").exists():
        run_command([tool_paths["bwa"], "index", "-p", host_reference_path, "-a", "bwtsw", host_reference_path], debug=args.debug)

    sai_path = temp_dir / "host_bwa.sai"
    aln_command = build_bwa_align_command(tool_paths, host_reference_path, sample_path, str(sai_path), args.threads)
    with sai_path.open("w", encoding="utf-8") as sai_handle:
        if args.debug:
            print("CMD:", " ".join(aln_command), ">", sai_path)
        completed = subprocess.run(aln_command, stdout=sai_handle, check=False)
        if completed.returncode != 0:
            raise SystemExit(f"Command failed with exit code {completed.returncode}: {' '.join(aln_command)}")

    samse_command = build_bwa_samse_command(tool_paths, host_reference_path, sai_path, sample_path, 1)
    with host_sam_path.open("w", encoding="utf-8") as sam_handle:
        if args.debug:
            print("CMD:", " ".join(samse_command), ">", host_sam_path)
        completed = subprocess.run(samse_command, stdout=sam_handle, check=False)
        if completed.returncode != 0:
            raise SystemExit(f"Command failed with exit code {completed.returncode}: {' '.join(samse_command)}")

    expand_xa_hits(host_sam_path)
    unmapped_num = generate_unmapped_reads_from_sam(host_sam_path, unmapped_path)
    mapped_num = max(count_fasta_sequences(sample_path) - unmapped_num, 0)
    return unmapped_path, mapped_num


def generate_aligned_contigs(sample_path: str, sam_path: str | Path, db_path: str, args, temp_dir: Path) -> tuple[Path, int]:
    tool_paths = tool_path_map()
    if "samtools" not in tool_paths:
        raise SystemExit("samtools is required for aligned-contig generation")

    reference_path = resolve_reference_path(args.reference, db_path)
    seq_info_path = resolve_seq_info_path(db_path)
    sample_name = Path(sample_path).name
    bam_path = temp_dir / f"{sample_name}.bam"
    sorted_bam_path = temp_dir / f"{sample_name}.sorted.bam"
    sorted_prefix = temp_dir / f"{sample_name}.sorted"
    pre_pileup_path = temp_dir / f"{sample_name}.pre.pileup"
    pileup_path = temp_dir / f"{sample_name}.pileup"
    aligned_path = temp_dir / f"{sample_name}.aligned.fa"

    view_command = [
        tool_paths["samtools"],
        "view",
        "-bt",
        f"{reference_path}.fai",
        str(sam_path),
    ]
    with bam_path.open("wb") as bam_handle:
        if args.debug:
            print("CMD:", " ".join(view_command), ">", bam_path)
        completed = subprocess.run(view_command, stdout=bam_handle, check=False)
        if completed.returncode != 0:
            raise SystemExit(f"Command failed with exit code {completed.returncode}: {' '.join(view_command)}")

    sort_command = [
        tool_paths["samtools"],
        "sort",
        str(bam_path),
        str(sorted_prefix),
    ]
    run_command(sort_command, debug=args.debug)

    if not sorted_bam_path.exists():
        raise SystemExit(f"samtools sort did not produce expected output: {sorted_bam_path}")

    mpileup_command = [
        tool_paths["samtools"],
        "mpileup",
        "-f",
        reference_path,
        str(sorted_bam_path),
    ]
    with pre_pileup_path.open("w", encoding="utf-8") as pileup_handle:
        if args.debug:
            print("CMD:", " ".join(mpileup_command), ">", pre_pileup_path)
        completed = subprocess.run(mpileup_command, stdout=pileup_handle, check=False)
        if completed.returncode != 0:
            raise SystemExit(f"Command failed with exit code {completed.returncode}: {' '.join(mpileup_command)}")

    filter_pileup_by_coverage(pre_pileup_path, seq_info_path, 0.3, pileup_path)
    contig_num = pileup_to_contig(reference_path, pileup_path, aligned_path, 40, 0, "ALIGNED")
    if contig_num == 0:
        aligned_path.write_text("", encoding="utf-8")

    return aligned_path, count_fasta_sequences(aligned_path)


def parse_kmer_values(kmer_range: str) -> list[int]:
    start_text, end_text = kmer_range.split("-", 1)
    start = int(start_text)
    end = int(end_text)
    if start > end:
        raise ValueError(f"Invalid k-mer range: {kmer_range}")
    values = [value for value in range(start, end + 1) if value % 2 == 1]
    if not values:
        raise ValueError(f"No odd k-mer values in range: {kmer_range}")
    return values


def summarize_contigs(contigs_path: Path) -> dict[str, float]:
    if not contigs_path.exists() or contigs_path.stat().st_size == 0:
        return {}

    lengths = []
    total_bases = 0
    total_ns = 0
    for record in read_fasta_records(contigs_path):
        length = len(record.sequence)
        if length == 0:
            continue
        lengths.append(length)
        total_bases += length
        total_ns += record.sequence.upper().count("N")

    if not lengths:
        return {}

    lengths.sort()
    cumulative = 0
    n50 = 0
    for length in lengths:
        cumulative += length
        if cumulative >= total_bases / 2:
            n50 = length
            break

    return {
        "numSeqs": float(len(lengths)),
        "numBases": float(total_bases),
        "numOK": float(total_bases - total_ns),
        "numNs": float(total_ns),
        "minLen": float(lengths[0]),
        "avgLen": float(total_bases) / len(lengths),
        "maxLen": float(lengths[-1]),
        "n50": float(n50),
    }


def prepare_assembly_input(sample_path: str, temp_dir: Path) -> Path:
    assembly_input_path = Path(sample_path)
    assembly_file_type = detect_file_type(sample_path)

    if assembly_input_path.suffix.lower() not in {".fa", ".fasta", ".fq", ".fastq", ".fa.gz", ".fasta.gz", ".fq.gz", ".fastq.gz"}:
        normalized_suffix = ".fq" if assembly_file_type == "fastq" else ".fa"
        normalized_input = temp_dir / f"{assembly_input_path.name}{normalized_suffix}"
        normalized_input.write_text(assembly_input_path.read_text(encoding="utf-8"), encoding="utf-8")
        assembly_input_path = normalized_input

    return assembly_input_path


def run_spades_assembly(sample_path: str, args, temp_dir: Path, data_type: str) -> tuple[Path, int]:
    tool_paths = tool_path_map()
    if "spades.py" not in tool_paths:
        raise SystemExit("spades.py is required for the Python de novo assembly stage")

    sample_name = Path(sample_path).name
    spades_dir = temp_dir / f"{sample_name}.spades"
    assembled_path = temp_dir / f"{sample_name}.assembled.fa"
    assembly_input_path = prepare_assembly_input(sample_path, temp_dir)

    if spades_dir.exists():
        shutil.rmtree(spades_dir, ignore_errors=True)

    command = [
        tool_paths["spades.py"],
        "--only-assembler",
        "-t",
        str(args.threads),
        "-o",
        str(spades_dir),
        "-s",
        str(assembly_input_path),
    ]
    if data_type == "mRNA":
        command.insert(1, "--rna")
    else:
        command.extend(["-k", ",".join(str(value) for value in parse_kmer_values(args.kmer_range))])

    run_command(command, debug=args.debug)

    contigs_path = spades_dir / "contigs.fasta"
    if not contigs_path.exists():
        assembled_path.write_text("", encoding="utf-8")
        return assembled_path, 0

    assembled_path.write_text(contigs_path.read_text(encoding="utf-8"), encoding="utf-8")
    return assembled_path, count_fasta_sequences(assembled_path)


def run_velvet(sample_path: Path, output_dir: Path, kmer_length: int, cov_cutoff: int, tool_paths: dict[str, str], debug: bool) -> Path | None:
    file_type = detect_file_type(sample_path)
    log_path = output_dir.parent / "velvet.log"
    shutil.rmtree(output_dir, ignore_errors=True)

    velveth_command = [
        tool_paths["velveth"],
        str(output_dir),
        str(kmer_length),
        f"-{file_type}",
        str(sample_path),
    ]
    if run_logged_command(velveth_command, log_path, debug=debug) != 0:
        shutil.rmtree(output_dir, ignore_errors=True)
        return None

    velvetg_command = [
        tool_paths["velvetg"],
        str(output_dir),
        "-cov_cutoff",
        str(cov_cutoff),
        "-min_contig_lgth",
        "30",
    ]
    if run_logged_command(velvetg_command, log_path, debug=debug) != 0:
        shutil.rmtree(output_dir, ignore_errors=True)
        return None

    contigs_path = output_dir / "contigs.fa"
    if not contigs_path.exists():
        shutil.rmtree(output_dir, ignore_errors=True)
        return None

    return contigs_path


def run_velvet_assembly(sample_path: str, args, temp_dir: Path, data_type: str) -> tuple[Path, int]:
    tool_paths = tool_path_map()
    if "velveth" not in tool_paths or "velvetg" not in tool_paths:
        raise SystemExit("velveth and velvetg are required for the Velvet de novo assembly stage")

    sample_name = Path(sample_path).name
    assembled_path = temp_dir / f"{sample_name}.assembled.fa"
    assembly_input_path = prepare_assembly_input(sample_path, temp_dir)
    velvet_dir = temp_dir / f"{sample_name}.velvet"
    velvet_input_path = assembly_input_path

    if args.rm_dup:
        deduplicated_input = temp_dir / f"{assembly_input_path.name}.uniq"
        remove_duplicate_reads(assembly_input_path, deduplicated_input)
        velvet_input_path = deduplicated_input

    if data_type == "mRNA":
        best_kmer = 31
        best_coverage = 10
    else:
        kmer_values = parse_kmer_values(args.kmer_range)
        best_kmer = kmer_values[0]
        best_coverage = 5
        best_objective = -1.0

        for kmer_value in kmer_values:
            candidate_dir = temp_dir / f"{sample_name}.velvet.{kmer_value}.5"
            contigs_path = run_velvet(
                velvet_input_path,
                candidate_dir,
                kmer_value,
                5,
                tool_paths,
                debug=args.debug,
            )
            stats = summarize_contigs(contigs_path) if contigs_path is not None else {}
            objective = stats.get("maxLen", -1.0)
            if objective > best_objective:
                best_objective = objective
                best_kmer = kmer_value
            shutil.rmtree(candidate_dir, ignore_errors=True)

        for coverage_value in range(7, 26):
            candidate_dir = temp_dir / f"{sample_name}.velvet.{best_kmer}.{coverage_value}"
            contigs_path = run_velvet(
                velvet_input_path,
                candidate_dir,
                best_kmer,
                coverage_value,
                tool_paths,
                debug=args.debug,
            )
            stats = summarize_contigs(contigs_path) if contigs_path is not None else {}
            objective = stats.get("maxLen", -1.0)
            if objective > best_objective:
                best_objective = objective
                best_coverage = coverage_value
            shutil.rmtree(candidate_dir, ignore_errors=True)

        if args.debug:
            print(
                "Velvet best parameters:",
                f"kmer={best_kmer}",
                f"coverage={best_coverage}",
                f"objective=maxLen",
            )

    final_contigs = run_velvet(
        velvet_input_path,
        velvet_dir,
        best_kmer,
        best_coverage,
        tool_paths,
        debug=args.debug,
    )
    if final_contigs is None:
        assembled_path.write_text("", encoding="utf-8")
        return assembled_path, 0

    assembled_path.write_text(final_contigs.read_text(encoding="utf-8"), encoding="utf-8")
    return assembled_path, count_fasta_sequences(assembled_path)


def run_denovo_assembly(sample_path: str, args, temp_dir: Path, data_type: str) -> tuple[Path, int]:
    if args.assembler == "spades":
        return run_spades_assembly(sample_path, args, temp_dir, data_type)
    if args.assembler == "velvet":
        return run_velvet_assembly(sample_path, args, temp_dir, data_type)
    raise SystemExit(f"Unsupported assembler: {args.assembler}")


def remove_redundant_contigs(input_path: str | Path, output_path: str | Path) -> tuple[int, int]:
    records = read_fasta_records(input_path)
    if not records:
        Path(output_path).write_text("", encoding="utf-8")
        return 0, 0

    indexed_records = [(index, record) for index, record in enumerate(records) if record.sequence]
    sorted_records = sorted(indexed_records, key=lambda item: (-len(item[1].sequence), item[0]))

    kept: list[tuple[int, FastaRecord]] = []
    kept_sequences: list[str] = []
    kept_reverse_complements: list[str] = []

    for index, record in sorted_records:
        sequence = record.sequence.upper()
        reverse_sequence = reverse_complement(sequence).upper()

        redundant = False
        for kept_sequence, kept_reverse_sequence in zip(kept_sequences, kept_reverse_complements):
            if (
                sequence in kept_sequence
                or sequence in kept_reverse_sequence
                or reverse_sequence in kept_sequence
                or reverse_sequence in kept_reverse_sequence
            ):
                redundant = True
                break

        if redundant:
            continue

        kept.append((index, record))
        kept_sequences.append(sequence)
        kept_reverse_complements.append(reverse_sequence)

    kept.sort(key=lambda item: item[0])
    output_records = [record for _index, record in kept]
    write_fasta_records(output_path, output_records)
    return len(records), len(output_records)
