from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from virusdetect import __db_version__, __version__
from virusdetect.db import (
    bundle_database_archive,
    install_database_archive,
    resolve_database_download_source,
    resolve_database_location,
    resolve_database_target_dir,
    verify_database_files,
)
from virusdetect.pipeline import run_pipeline
from virusdetect.tools import check_tools, install_hint_text, missing_required_tools


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="virusdetect", description="VirusDetect v2 command line interface")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(dest="command", required=True)

    build_run_parser(subparsers)
    build_db_parser(subparsers)
    build_tools_parser(subparsers)

    return parser


def build_run_parser(subparsers) -> None:
    parser = subparsers.add_parser("run", help="Run the Python pipeline scaffold")
    parser.add_argument("input_paths", nargs="+", help="Input read or contig file(s)")
    parser.add_argument("--backend", choices=["legacy", "python"], default="python", help="Execution backend")
    parser.add_argument("--db-path", help="Path to a VirusDetect database directory")
    parser.add_argument("--reference", default="vrl_plant", help="Reference database name")
    parser.add_argument("--host-reference", help="Host reference database path or name")
    parser.add_argument("-o", "--output", help="Output directory")
    parser.add_argument("-t", "--threads", type=int, default=8, help="Number of threads")
    parser.add_argument("--read-length", help="sRNA length selector, for example 21-23:25:27:32-34")
    parser.add_argument("--kmer-range", default="9-23", help="De novo assembly k-mer range")
    parser.add_argument("--hisat-dist", type=int, default=5, help="Maximum edit distance surrogate for HISAT2 scoring")
    parser.add_argument("--word-size", type=int, default=11, help="Minimum word size for nucleotide search")
    parser.add_argument("--exp-value", type=float, default=1e-5, help="Maximum e-value for nucleotide search")
    parser.add_argument("--exp-valuex", type=float, help="Override translated search e-value")
    parser.add_argument("--percent-identity", type=float, default=25, help="Minimum translated hit identity percentage")
    parser.add_argument("--hsp-cover", type=float, default=0.75, help="Minimum HSP coverage")
    parser.add_argument("--coverage-cutoff", type=float, default=0.1, help="Reference coverage cutoff")
    parser.add_argument("--depth-cutoff", type=float, default=5.0, help="Reference depth cutoff")
    parser.add_argument("--norm-depth-cutoff", type=float, default=5.0, help="Normalized reference depth cutoff")
    parser.add_argument("--diff-ratio", type=float, default=0.25, help="Maximum unique-contig ratio for redundant references")
    parser.add_argument(
        "--diff-contig-cover",
        type=float,
        default=0.5,
        help="Maximum unique reference coverage ratio for redundant references",
    )
    parser.add_argument(
        "--diff-contig-length",
        type=int,
        default=100,
        help="Maximum unique covered bases for redundant references",
    )
    parser.add_argument("--siRNA-percent", type=float, default=0.5, help="21-22 nt siRNA enrichment cutoff")
    parser.add_argument("--check-only", action="store_true", help="Validate runtime inputs without running analysis")
    parser.add_argument("--prepare-only", action="store_true", help="Prepare sample inputs and emit analysis plans")
    parser.add_argument(
        "--stop-after",
        choices=[
            "prepare_sample",
            "align_to_virus_reference",
            "host_subtraction",
            "de_novo_assembly",
            "generate_aligned_contigs",
            "combine_contigs",
            "remove_redundant_contigs",
            "identify_virus",
        ],
        help="Stop the Python backend after a specific stage",
    )
    parser.add_argument("--keep-temp", action="store_true", help="Keep temporary working files")
    parser.add_argument("--debug", action="store_true", help="Enable verbose backend logging")
    parser.set_defaults(handler=handle_run)


def build_db_parser(subparsers) -> None:
    parser = subparsers.add_parser("db", help="Database management commands")
    db_subparsers = parser.add_subparsers(dest="db_command", required=True)

    path_parser = db_subparsers.add_parser("path", help="Print the resolved database path")
    path_parser.add_argument("--target", action="store_true", help="Print the default install target path")
    path_parser.set_defaults(handler=handle_db)

    verify_parser = db_subparsers.add_parser("verify", help="Verify required database files")
    verify_parser.add_argument("--path", dest="db_path", help="Database path to verify")
    verify_parser.set_defaults(handler=handle_db)

    download_parser = db_subparsers.add_parser("download", help="Download and install a database archive")
    download_parser.add_argument("--path", dest="db_path", help="Install database into this directory")
    download_parser.add_argument("--url", help="Database archive URL")
    download_parser.add_argument("--sha256", help="Expected SHA256 for the archive")
    download_parser.add_argument("--sha256-url", help="URL pointing to a SHA256 checksum file")
    download_parser.add_argument("--repo", help="GitHub repository that hosts database release assets")
    download_parser.add_argument("--release-tag", help="Git tag of the release that hosts database assets")
    download_parser.add_argument("--asset-name", help="Database archive asset name")
    download_parser.add_argument("--sha256-asset-name", help="Checksum asset name")
    download_parser.add_argument("--db-version", default=__db_version__, help="Default database bundle version label")
    download_parser.set_defaults(handler=handle_db)

    bundle_parser = db_subparsers.add_parser("bundle", help="Create a distributable database archive")
    bundle_parser.add_argument("--path", dest="db_path", help="Source database path to package")
    bundle_parser.add_argument("--output", help="Output `.tar.gz` path")
    bundle_parser.add_argument("--db-version", default=__version__, help="Database bundle version label")
    bundle_parser.add_argument("--name", default="virusdetect", help="Database bundle name")
    bundle_parser.set_defaults(handler=handle_db)


def build_tools_parser(subparsers) -> None:
    parser = subparsers.add_parser("tools", help="External tool checks")
    tools_subparsers = parser.add_subparsers(dest="tools_command", required=True)

    check_parser = tools_subparsers.add_parser("check", help="Check whether external tools are available")
    check_parser.add_argument("--json", action="store_true", help="Emit machine-readable JSON")
    check_parser.add_argument("--strict", action="store_true", help="Exit non-zero when required tools are missing")
    check_parser.set_defaults(handler=handle_tools)

    install_parser = tools_subparsers.add_parser("install-hint", help="Print recommended pixi/conda install commands")
    install_parser.set_defaults(handler=handle_tools_install_hint)


def handle_run(args) -> int:
    return run_pipeline(args)


def handle_db(args) -> int:
    if args.db_command == "path":
        if args.target:
            print(resolve_database_target_dir())
            return 0

        location = resolve_database_location()
        if location is None:
            print("Database not found")
            return 1

        print(location.path)
        return 0

    if args.db_command == "verify":
        if args.db_path:
            db_path = args.db_path
        else:
            location = resolve_database_location()
            db_path = location.path if location is not None else resolve_database_target_dir()

        missing = verify_database_files(db_path)
        if missing:
            print(f"Database verification failed: {db_path}")
            for section, file_names in sorted(missing.items()):
                print(f"{section}: {', '.join(file_names)}")
            return 1

        print(f"Database verification passed: {db_path}")
        return 0

    if args.db_command == "download":
        target_dir = resolve_database_target_dir(args.db_path)
        source = resolve_database_download_source(
            url=args.url,
            sha256=args.sha256,
            sha256_url=args.sha256_url,
            release_repo=args.repo,
            release_tag=args.release_tag,
            asset_name=args.asset_name,
            sha256_asset_name=args.sha256_asset_name,
            db_version=args.db_version,
        )
        installed_path = install_database_archive(
            url=source.url,
            destination_dir=target_dir,
            sha256=source.sha256,
            sha256_url=source.sha256_url,
        )
        print(f"Database source: {source.url}")
        print(f"Database installed to: {installed_path}")
        return 0

    if args.db_command == "bundle":
        if args.db_path:
            source_path = args.db_path
        else:
            location = resolve_database_location()
            if location is None:
                raise SystemExit("No VirusDetect database was found to bundle")
            source_path = location.path

        output_path = args.output or str(Path("dist") / f"virusdetect-db-{args.db_version}.tar.gz")
        bundle = bundle_database_archive(
            source_dir=source_path,
            destination_archive=output_path,
            db_version=args.db_version,
            db_name=args.name,
        )
        print(f"Database bundle created: {bundle.archive_path}")
        print(f"SHA256 file written: {bundle.sha256_path}")
        return 0

    raise SystemExit(f"Unknown db command: {args.db_command}")


def handle_tools(args) -> int:
    statuses = check_tools()
    if args.json:
        print(
            json.dumps(
                [
                    {
                        "name": status.name,
                        "required": status.required,
                        "found": status.found,
                        "path": status.path,
                        "source": status.source,
                    }
                    for status in statuses
                ],
                indent=2,
            )
        )
    else:
        for status in statuses:
            label = "required" if status.required else "optional"
            state = "OK" if status.found else "MISSING"
            path = status.path or "-"
            source = status.source or "-"
            print(f"{status.name}\t{label}\t{state}\t{source}\t{path}")

    missing = missing_required_tools(statuses)
    return 1 if args.strict and missing else 0


def handle_tools_install_hint(_args) -> int:
    print(install_hint_text())
    return 0


def normalize_argv(argv: list[str]) -> list[str]:
    if argv and argv[0] == "--":
        return argv[1:]
    return argv


def main(argv: list[str] | None = None) -> None:
    normalized_argv = normalize_argv(sys.argv[1:] if argv is None else argv)
    parser = build_parser()
    args = parser.parse_args(normalized_argv)
    raise SystemExit(args.handler(args))
