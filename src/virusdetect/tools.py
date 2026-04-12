CORE_CONDA_PACKAGES = ("bwa", "samtools", "blast")
OPTIONAL_CONDA_PACKAGES = ("hisat2", "spades", "velvet")


def required_runtime_tools(*, data_type: str, host_reference: str | None, assembler: str) -> tuple[str, ...]:
    required = ["bwa", "samtools", "blastn", "blastx", "makeblastdb"]

    if host_reference and data_type == "mRNA":
        required.extend(["hisat2", "hisat2-build"])

    if assembler == "spades":
        required.append("spades.py")
    elif assembler == "velvet":
        required.extend(["velveth", "velvetg"])
    else:
        raise ValueError(f"Unsupported assembler: {assembler}")

    return tuple(required)

from virusdetect.runtime import missing_required_tools, resolve_tools


def check_tools():
    return resolve_tools()


def install_hint_text() -> str:
    core_packages = " ".join(CORE_CONDA_PACKAGES)
    optional_packages = " ".join(OPTIONAL_CONDA_PACKAGES)
    lines = [
        "Recommended pixi setup for the Python backend:",
        f"  pixi add {core_packages} {optional_packages}",
        "",
        "Equivalent conda / mamba setup:",
        f"  mamba install -c conda-forge -c bioconda {core_packages} {optional_packages}",
        "",
        "Notes:",
        "  - `bwa`, `samtools`, and `blast` are required for the current Python v2 stages.",
        "  - `hisat2` is used for mRNA host subtraction when `--host-reference` is enabled.",
        "  - `spades` is the default de novo assembler for the current Python pipeline.",
        "  - `velvet` is available as an alternate assembler via `--assembler velvet` for v1-style comparison.",
        "  - the default `main` workflow no longer needs Perl packages.",
        "  - for the historical Perl workflow, use the `v1` branch and install its dependencies there.",
    ]
    return "\n".join(lines)
