from __future__ import annotations

import os
import shlex
import subprocess
from pathlib import Path

from virusdetect.db import resolve_reference_path
from virusdetect.runtime import repo_root, resolve_tool

LEGACY_PERL_MODULES = (
    "Bio::SeqIO",
    "Bio::Graphics",
)


def legacy_script_path() -> Path:
    return repo_root() / "virus_detect.pl"


def legacy_identifier_path() -> Path:
    return repo_root() / "bin" / "virus_identify.pl"


def build_legacy_command(args, input_paths: list[Path]) -> list[str]:
    command = [
        "perl",
        str(legacy_script_path()),
        "--reference",
        args.reference,
        "--thread_num",
        str(args.threads),
        "--kmer_range",
        args.kmer_range,
        "--word_size",
        str(args.word_size),
        "--exp_value",
        str(args.exp_value),
        "--percent_identity",
        str(args.percent_identity),
        "--hsp_cover",
        str(args.hsp_cover),
        "--coverage_cutoff",
        str(args.coverage_cutoff),
        "--depth_cutoff",
        str(args.depth_cutoff),
        "--norm_depth_cutoff",
        str(args.norm_depth_cutoff),
        "--diff-ratio",
        str(args.diff_ratio),
        "--diff-contig-cover",
        str(args.diff_contig_cover),
        "--diff-contig-length",
        str(args.diff_contig_length),
        "--siRNA_percent",
        str(args.siRNA_percent),
    ]

    if args.host_reference:
        command.extend(["--host_reference", args.host_reference])
    if args.read_length:
        command.extend(["--read_length", args.read_length])
    if args.exp_valuex is not None:
        command.extend(["--exp-valuex", str(args.exp_valuex)])
    if args.debug:
        command.append("--debug")

    command.extend(str(path.resolve()) for path in input_paths)
    return command


def run_legacy_pipeline(args, input_paths: list[Path], workdir: Path) -> int:
    workdir.mkdir(parents=True, exist_ok=True)
    command = build_legacy_command(args, input_paths)
    env = legacy_execution_env(workdir)

    if args.debug:
        print("Legacy backend command:")
        print(" ".join(command))

    completed = subprocess.run(command, cwd=str(workdir), check=False, env=env)
    return completed.returncode


def build_legacy_identifier_command(args, db_path: str, sample_path: Path, contig_path: Path) -> list[str]:
    reference_path = resolve_reference_path(args.reference, db_path)
    command = [
        "perl",
        str(legacy_identifier_path()),
        "--reference",
        reference_path,
        "--word-size",
        str(args.word_size),
        "--exp-value",
        str(args.exp_value),
        "--cpu-num",
        str(args.threads),
        "--percent-identity",
        str(args.percent_identity),
        "--hsp-cover",
        str(args.hsp_cover),
        "--coverage-cutoff",
        str(args.coverage_cutoff),
        "--depth-cutoff",
        str(args.depth_cutoff),
        "--norm-depth-cutoff",
        str(args.norm_depth_cutoff),
        "--diff-ratio",
        str(args.diff_ratio),
        "--diff-contig-cover",
        str(args.diff_contig_cover),
        "--diff-contig-length",
        str(args.diff_contig_length),
        "--siRNA-percent",
        str(args.siRNA_percent),
        str(sample_path.resolve()),
        str(contig_path.resolve()),
    ]
    if args.exp_valuex is not None:
        command[10:10] = ["--exp-valuex", str(args.exp_valuex)]
    if args.debug:
        command.insert(-2, "--debug")
    return command


def run_legacy_identifier(args, db_path: str, sample_path: Path, contig_path: Path, workdir: Path) -> int:
    workdir.mkdir(parents=True, exist_ok=True)
    command = build_legacy_identifier_command(args, db_path, sample_path, contig_path)
    env = legacy_execution_env(workdir)
    if args.debug:
        print("Legacy identifier command:")
        print(" ".join(command))
    completed = subprocess.run(command, cwd=str(workdir), check=False, env=env)
    return completed.returncode


def missing_legacy_perl_modules() -> list[str]:
    missing: list[str] = []
    perl_lib = repo_root() / "bin"
    for module_name in LEGACY_PERL_MODULES:
        completed = subprocess.run(
            ["perl", "-I", str(perl_lib), f"-M{module_name}", "-e", "1"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
        if completed.returncode != 0:
            missing.append(module_name)
    return missing


def legacy_execution_env(workdir: Path) -> dict[str, str]:
    env = os.environ.copy()
    bridge_dir = create_legacy_tool_bridge(workdir)
    env["VIRUSDETECT_TOOL_DIR"] = str(bridge_dir)
    env["PATH"] = str(bridge_dir) + os.pathsep + env.get("PATH", "")
    return env


def create_legacy_tool_bridge(workdir: Path) -> Path:
    bridge_dir = workdir / ".virusdetect-legacy-tools"
    bridge_dir.mkdir(parents=True, exist_ok=True)

    resolved_tools = _resolve_legacy_tool_paths()
    _write_proxy_script(bridge_dir / "bwa", resolved_tools["bwa"])
    _write_proxy_script(bridge_dir / "samtools", resolved_tools["samtools"])
    _write_megablast_shim(bridge_dir / "megablast", resolved_tools["blastn"])
    _write_blastall_shim(bridge_dir / "blastall", resolved_tools["blastx"])
    _write_formatdb_shim(bridge_dir / "formatdb", resolved_tools["makeblastdb"])
    return bridge_dir


def _resolve_legacy_tool_paths() -> dict[str, str]:
    resolved: dict[str, str] = {}
    for tool_name in ("bwa", "samtools", "blastn", "blastx", "makeblastdb"):
        path, _source = resolve_tool(tool_name)
        if path is None:
            raise SystemExit(f"Required tool for legacy bridge was not found: {tool_name}")
        resolved[tool_name] = path
    return resolved


def _write_proxy_script(path: Path, target: str) -> None:
    path.write_text(
        "#!/usr/bin/env bash\n"
        f"exec {shlex.quote(target)} \"$@\"\n",
        encoding="utf-8",
    )
    path.chmod(0o755)


def _write_megablast_shim(path: Path, blastn_path: str) -> None:
    path.write_text(
        "#!/usr/bin/env python3\n"
        "import os\n"
        "import sys\n"
        f"BLASTN = {blastn_path!r}\n"
        "translated = [BLASTN, '-task', 'megablast']\n"
        "gapopen = None\n"
        "gapextend = None\n"
        "args = sys.argv[1:]\n"
        "index = 0\n"
        "while index < len(args):\n"
        "    arg = args[index]\n"
        "    if arg == '-i' and index + 1 < len(args):\n"
        "        translated.extend(['-query', args[index + 1]])\n"
        "        index += 2\n"
        "    elif arg == '-d' and index + 1 < len(args):\n"
        "        translated.extend(['-db', args[index + 1]])\n"
        "        index += 2\n"
        "    elif arg == '-o' and index + 1 < len(args):\n"
        "        translated.extend(['-out', args[index + 1]])\n"
        "        index += 2\n"
        "    elif arg == '-F' and index + 1 < len(args):\n"
        "        translated.extend(['-dust', 'no' if args[index + 1] == 'F' else 'yes'])\n"
        "        index += 2\n"
        "    elif arg == '-a' and index + 1 < len(args):\n"
        "        translated.extend(['-num_threads', args[index + 1]])\n"
        "        index += 2\n"
        "    elif arg == '-W' and index + 1 < len(args):\n"
        "        translated.extend(['-word_size', args[index + 1]])\n"
        "        index += 2\n"
        "    elif arg == '-q' and index + 1 < len(args):\n"
        "        translated.extend(['-penalty', args[index + 1]])\n"
        "        index += 2\n"
        "    elif arg == '-G' and index + 1 < len(args):\n"
        "        gapopen = args[index + 1]\n"
        "        index += 2\n"
        "    elif arg == '-E' and index + 1 < len(args):\n"
        "        gapextend = args[index + 1]\n"
        "        index += 2\n"
        "    elif arg == '-b' and index + 1 < len(args):\n"
        "        translated.extend(['-max_target_seqs', args[index + 1]])\n"
        "        index += 2\n"
        "    elif arg == '-e' and index + 1 < len(args):\n"
        "        translated.extend(['-evalue', args[index + 1]])\n"
        "        index += 2\n"
        "    elif arg == '-S' and index + 1 < len(args):\n"
        "        strand = {'1': 'plus', '2': 'minus', '3': 'both'}.get(args[index + 1])\n"
        "        if strand is not None:\n"
        "            translated.extend(['-strand', strand])\n"
        "        index += 2\n"
        "    else:\n"
        "        translated.append(arg)\n"
        "        index += 1\n"
        "if gapopen is not None and gapextend is not None:\n"
        "    try:\n"
        "        if int(gapopen) >= 0 and int(gapextend) >= 0:\n"
        "            translated.extend(['-gapopen', gapopen, '-gapextend', gapextend])\n"
        "    except ValueError:\n"
        "        pass\n"
        "os.execv(translated[0], translated)\n",
        encoding="utf-8",
    )
    path.chmod(0o755)


def _write_blastall_shim(path: Path, blastx_path: str) -> None:
    path.write_text(
        "#!/usr/bin/env python3\n"
        "import os\n"
        "import sys\n"
        f"BLASTX = {blastx_path!r}\n"
        "program = 'blastx'\n"
        "translated = [BLASTX]\n"
        "args = sys.argv[1:]\n"
        "index = 0\n"
        "while index < len(args):\n"
        "    arg = args[index]\n"
        "    if arg == '-p' and index + 1 < len(args):\n"
        "        program = args[index + 1]\n"
        "        index += 2\n"
        "    elif arg == '-i' and index + 1 < len(args):\n"
        "        translated.extend(['-query', args[index + 1]])\n"
        "        index += 2\n"
        "    elif arg == '-d' and index + 1 < len(args):\n"
        "        translated.extend(['-db', args[index + 1]])\n"
        "        index += 2\n"
        "    elif arg == '-o' and index + 1 < len(args):\n"
        "        translated.extend(['-out', args[index + 1]])\n"
        "        index += 2\n"
        "    elif arg == '-F' and index + 1 < len(args):\n"
        "        translated.extend(['-seg', 'no' if args[index + 1] == 'F' else 'yes'])\n"
        "        index += 2\n"
        "    elif arg == '-a' and index + 1 < len(args):\n"
        "        translated.extend(['-num_threads', args[index + 1]])\n"
        "        index += 2\n"
        "    elif arg == '-e' and index + 1 < len(args):\n"
        "        translated.extend(['-evalue', args[index + 1]])\n"
        "        index += 2\n"
        "    else:\n"
        "        translated.append(arg)\n"
        "        index += 1\n"
        "if program != 'blastx':\n"
        "    sys.stderr.write(f'blastall shim only supports -p blastx, got {program}\\n')\n"
        "    sys.exit(2)\n"
        "os.execv(translated[0], translated)\n",
        encoding="utf-8",
    )
    path.chmod(0o755)


def _write_formatdb_shim(path: Path, makeblastdb_path: str) -> None:
    path.write_text(
        "#!/usr/bin/env python3\n"
        "import os\n"
        "import sys\n"
        f"MAKEBLASTDB = {makeblastdb_path!r}\n"
        "input_path = None\n"
        "dbtype = 'nucl'\n"
        "args = sys.argv[1:]\n"
        "index = 0\n"
        "while index < len(args):\n"
        "    arg = args[index]\n"
        "    if arg == '-i' and index + 1 < len(args):\n"
        "        input_path = args[index + 1]\n"
        "        index += 2\n"
        "    elif arg == '-p' and index + 1 < len(args):\n"
        "        dbtype = 'prot' if args[index + 1] == 'T' else 'nucl'\n"
        "        index += 2\n"
        "    else:\n"
        "        index += 1\n"
        "if input_path is None:\n"
        "    sys.stderr.write('formatdb shim requires -i <input>\\n')\n"
        "    sys.exit(2)\n"
        "os.execv(MAKEBLASTDB, [MAKEBLASTDB, '-in', input_path, '-dbtype', dbtype])\n",
        encoding="utf-8",
    )
    path.chmod(0o755)
