from __future__ import annotations

import json
import shutil
from dataclasses import asdict, dataclass
from pathlib import Path

from virusdetect.db import resolve_database_location, verify_database_files
from virusdetect.fasta import read_fasta_records, write_fasta_records
from virusdetect.identify import run_python_identifier
from virusdetect.sample import count_sequences, detect_data_type, detect_file_type, filter_reads_by_length
from virusdetect.stages import (
    generate_aligned_contigs,
    remove_redundant_contigs,
    run_denovo_assembly,
    run_host_subtraction,
    run_virus_alignment,
)
from virusdetect.tools import check_tools, missing_required_tools


@dataclass(frozen=True)
class SamplePlan:
    input_path: str
    sample_name: str
    file_type: str
    data_type: str
    sequence_count: int
    temp_dir: str
    prepared_input: str
    read_length: str | None
    stages: list[str]


def resolve_sample_paths(input_path: Path, output_dir: str | None) -> tuple[Path, Path]:
    sample_name = input_path.name
    if output_dir:
        sample_output_dir = Path(output_dir).resolve()
    else:
        sample_output_dir = Path.cwd()
    temp_dir = sample_output_dir / f"{sample_name}_temp"
    return sample_output_dir, temp_dir


def create_sample_plan(input_path: Path, args) -> SamplePlan:
    file_type = detect_file_type(input_path)
    data_type = detect_data_type(input_path)
    sequence_count = count_sequences(input_path)
    sample_output_dir, temp_dir = resolve_sample_paths(input_path, args.output)
    prepared_input = temp_dir / input_path.name
    effective_exp_valuex = args.exp_valuex if args.exp_valuex is not None else (1e-5 if data_type == "mRNA" else 1e-2)

    stages = [
        "detect_input",
        "prepare_sample",
        "align_to_virus_reference",
        "generate_aligned_contigs",
        "de_novo_assembly",
        "remove_redundant_contigs",
        "identify_virus",
    ]
    if args.host_reference:
        stages.insert(4, "host_subtraction")

    return SamplePlan(
        input_path=str(input_path.resolve()),
        sample_name=input_path.name,
        file_type=file_type,
        data_type=data_type,
        sequence_count=sequence_count,
        temp_dir=str(temp_dir),
        prepared_input=str(prepared_input),
        read_length=args.read_length,
        stages=stages + [f"blastx_evalue={effective_exp_valuex}"],
    )


def prepare_sample_input(sample_plan: SamplePlan, args) -> None:
    temp_dir = Path(sample_plan.temp_dir)
    temp_dir.mkdir(parents=True, exist_ok=True)
    prepared_input = Path(sample_plan.prepared_input)
    source_input = Path(sample_plan.input_path)

    if args.read_length:
        kept = filter_reads_by_length(args.read_length, source_input, prepared_input)
        if kept == 0:
            raise SystemExit(f"No reads remained after read-length filtering for {source_input}")
        return

    if prepared_input.exists():
        return
    shutil.copy2(source_input, prepared_input)


def write_sample_plan(
    sample_plan: SamplePlan,
    sample_output_dir: Path,
    db_path: str,
    missing_tools: list[str],
) -> Path:
    sample_output_dir.mkdir(parents=True, exist_ok=True)
    plan_path = sample_output_dir / f"{sample_plan.sample_name}.analysis_plan.json"
    plan_payload = asdict(sample_plan)
    plan_payload["database"] = db_path
    plan_payload["missing_required_tools"] = missing_tools
    with plan_path.open("w", encoding="utf-8") as handle:
        json.dump(plan_payload, handle, indent=2, sort_keys=True)
    return plan_path


def combine_contigs(aligned_path: Path, assembled_path: Path, combined_path: Path) -> int:
    combined_records = []
    for candidate in (aligned_path, assembled_path):
        if candidate.exists() and candidate.stat().st_size > 0:
            combined_records.extend(read_fasta_records(candidate))

    write_fasta_records(combined_path, combined_records)
    return len(combined_records)


def run_pipeline(args) -> int:
    input_paths = [Path(path) for path in args.input_paths]
    missing_inputs = [str(path) for path in input_paths if not path.exists()]
    if missing_inputs:
        raise SystemExit(f"Input file(s) not found: {', '.join(missing_inputs)}")

    if args.db_path:
        db_path = args.db_path
    else:
        location = resolve_database_location()
        if location is None:
            raise SystemExit(
                "No VirusDetect database was found. Run `virusdetect db verify` or `virusdetect db download` first."
            )
        db_path = location.path

    missing_db = verify_database_files(db_path)
    if missing_db:
        raise SystemExit(f"Database verification failed for {db_path}: {missing_db}")

    tool_statuses = check_tools()
    missing_tools = missing_required_tools(tool_statuses)
    overall_exit_code = 0 if not missing_tools else 1

    for input_path in input_paths:
        sample_plan = create_sample_plan(input_path, args)
        sample_output_dir, _temp_dir = resolve_sample_paths(input_path, args.output)
        plan_path = write_sample_plan(sample_plan, sample_output_dir, db_path, missing_tools)

        print(f"Input: {input_path}")
        print(f"Database: {db_path}")
        print(f"Threads: {args.threads}")
        print(f"File type: {sample_plan.file_type}")
        print(f"Data type: {sample_plan.data_type}")
        print(f"Sequences: {sample_plan.sequence_count}")
        print(f"Analysis plan: {plan_path}")
        if args.output:
            print(f"Output: {sample_output_dir}")

        if missing_tools:
            print("Missing required tools: " + ", ".join(missing_tools))
        else:
            print("Required tools: OK")

        if args.check_only:
            continue

        prepare_sample_input(sample_plan, args)
        print(f"Prepared input: {sample_plan.prepared_input}")

        if args.prepare_only and not args.keep_temp:
            print("Temporary files retained because later stages are not ported yet.")

        if not args.prepare_only:
            sam_path, mapped_num = run_virus_alignment(sample_plan.prepared_input, db_path, args, Path(sample_plan.temp_dir))
            print(f"Virus SAM: {sam_path}")
            print(f"Mapped reads after filtering: {mapped_num}")
            if args.stop_after == "align_to_virus_reference":
                continue
            aligned_path, aligned_contigs = generate_aligned_contigs(
                sample_plan.prepared_input,
                sam_path,
                db_path,
                args,
                Path(sample_plan.temp_dir),
            )
            print(f"Aligned contigs: {aligned_path}")
            print(f"Aligned contig count: {aligned_contigs}")
            if args.stop_after == "generate_aligned_contigs":
                continue
            assembly_input_path, host_mapped_num = run_host_subtraction(
                sample_plan.prepared_input,
                db_path,
                args,
                Path(sample_plan.temp_dir),
                sample_plan.file_type,
                sample_plan.data_type,
            )
            print(f"Assembly input: {assembly_input_path}")
            if args.host_reference:
                print(f"Host-mapped reads: {host_mapped_num}")
            if args.stop_after == "host_subtraction":
                continue
            assembled_path, assembled_contigs = run_denovo_assembly(
                str(assembly_input_path),
                args,
                Path(sample_plan.temp_dir),
                sample_plan.data_type,
            )
            print(f"Assembled contigs: {assembled_path}")
            print(f"Assembled contig count: {assembled_contigs}")
            if args.stop_after == "de_novo_assembly":
                continue
            combined_path = Path(sample_plan.temp_dir) / f"{Path(sample_plan.prepared_input).name}.combined.fa"
            combined_contigs = combine_contigs(aligned_path, assembled_path, combined_path)
            print(f"Combined contigs: {combined_path}")
            print(f"Combined contig count: {combined_contigs}")
            if args.stop_after == "combine_contigs":
                continue
            deduplicated_path = Path(sample_plan.temp_dir) / f"{Path(sample_plan.prepared_input).name}.combined.nr.fa"
            _raw_combined_count, deduplicated_contigs = remove_redundant_contigs(combined_path, deduplicated_path)
            print(f"Nonredundant contigs: {deduplicated_path}")
            print(f"Nonredundant contig count: {deduplicated_contigs}")
            if args.stop_after == "remove_redundant_contigs":
                continue
            if deduplicated_contigs == 0:
                print("No combined contigs were generated; skipping identification.")
                if args.stop_after == "identify_virus":
                    continue
            else:
                identify_result = run_python_identifier(
                    args,
                    db_path,
                    Path(sample_plan.prepared_input),
                    deduplicated_path,
                    sample_output_dir,
                    sample_plan.data_type,
                )
                print(f"Result directory: {identify_result.result_dir}")
                print(f"Known contigs: {identify_result.known_contig_count}")
                print(f"Novel contigs: {identify_result.novel_contig_count}")
                print(f"Undetermined contigs: {identify_result.undetermined_contig_count}")
                print(f"Summary TSV: {identify_result.summary_tsv}")
                if args.stop_after == "identify_virus":
                    continue

    if args.check_only or args.prepare_only or args.stop_after in {
        "align_to_virus_reference",
        "generate_aligned_contigs",
        "host_subtraction",
        "de_novo_assembly",
        "combine_contigs",
        "remove_redundant_contigs",
        "identify_virus",
    }:
        return overall_exit_code
    return overall_exit_code
