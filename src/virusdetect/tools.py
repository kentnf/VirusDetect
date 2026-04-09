CORE_CONDA_PACKAGES = ("bwa", "samtools", "blast")
OPTIONAL_CONDA_PACKAGES = ("hisat2", "spades")
LEGACY_CONDA_PACKAGES = ("perl", "perl-bioperl")

from virusdetect.runtime import missing_required_tools, resolve_tools


def check_tools():
    return resolve_tools()


def install_hint_text() -> str:
    core_packages = " ".join(CORE_CONDA_PACKAGES)
    optional_packages = " ".join(OPTIONAL_CONDA_PACKAGES)
    legacy_packages = " ".join(LEGACY_CONDA_PACKAGES)
    return "\n".join(
        [
            "Recommended pixi setup:",
            f"  pixi add {core_packages} {optional_packages} {legacy_packages}",
            "",
            "Equivalent conda / mamba setup:",
            f"  mamba install -c conda-forge -c bioconda {core_packages} {optional_packages} {legacy_packages}",
            "",
            "Notes:",
            "  - `bwa`, `samtools`, and `blast` are required for the current Python v2 stages.",
            "  - `hisat2` and `spades` are used by the current Python pipeline for host subtraction and de novo assembly.",
            "  - legacy Perl execution still needs BioPerl-related modules; `perl-bioperl` is the main Bioconda package.",
            "  - `Bio::Graphics` may still need extra platform-specific setup beyond the packages above.",
        ]
    )
