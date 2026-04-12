# VirusDetect v2

Current package version: `2.0.0a0`

This branch is the Python rewrite line for VirusDetect.

Legacy Perl releases remain available on the `v1` branch and the historical tags (`v1.2` through `v1.8.1`). The old Perl entrypoints have been removed from `main`; new development should target the Python package under `src/virusdetect`.

## Status

The v2 branch currently provides:

- a Python package skeleton
- a `virusdetect` command-line entrypoint
- database path, verification, and archive download commands
- external tool discovery checks
- a `run` command that can validate inputs, classify samples, and prepare sample working files
- Python implementations of sample preparation, virus-reference alignment, host subtraction, de novo assembly, aligned contig generation, combined-contig deduplication, and virus identification
- Python-native `blastn`/`blastx` identification outputs, including Perl-style contig-classification FASTA files, reference FASTA exports, summary HTML pages, undetermined-contig reports, and compatibility `SAM` / `XLS` outputs
- a Python-first `run` command for the v2 rewrite line

The default backend on `main` is `python`. The old Perl workflow is no longer runnable through the `virusdetect` CLI on `main`; use the `v1` branch for supported legacy execution.

Branch roles:

- `main`: active Python v2 development line
- `v1`: legacy Perl maintenance line
- `master`: frozen historical alias of the old v1 state; kept temporarily for compatibility

Historical standalone utilities under `tools/` are kept on `main` for separate upstream/downstream use, but they are not part of the supported Python runtime path.

This is still an alpha release line: the Python pipeline is functional through identification and reporting, but report parity and package-manager distribution are not finished yet.

## Migration Snapshot

Implemented in the Python v2 line:

- CLI entrypoint and staged pipeline execution
- database path resolution, verification, and archive download plumbing
- sample preparation and input validation
- virus-reference alignment
- host subtraction
- de novo assembly
- aligned contig generation
- combined-contig redundancy reduction
- Python-native `blastn` and `blastx` identification
- redundancy pruning and coverage/depth filtering similar to the legacy Perl logic
- report artifacts:
  - Perl-style release output set for `identify_virus`
  - compatibility `blastn.sam`, `blastx.sam`, `blastn.xls`, and `blastx.xls`
  - reference FASTA export
  - summary HTML with compact coverage maps and stable per-reference row anchors
  - undetermined-contig HTML reports with grouped headers, browser-native tooltips, candidate highlighting, stable anchors, `undetermined.html` as the no-match view, and `undetermined_blast.html` as the match-only view for the current sRNA path

Not yet ported from the legacy Perl workflow:

- exact `Bio::Graphics` parity and detailed historical graphical layout behavior
- full parity for the historical siRNA tables and highlighting presentation
- historical per-reference detail HTML pages are not part of the default v2 release output
- full parity for every historical HTML artifact and formatting detail
- final publication flow for package-manager distribution and future database refreshes

## Quick Start

### Recommended: pixi

Install the local development environment:

```bash
pixi install
```

The `pixi` environment now includes the core command-line tools needed for the current migration stages:

- `bwa`
- `samtools`
- `blast`
- `hisat2`
- `spades`
- `velvet`

Legacy-only extras are no longer installed on `main`. If you still need the historical Perl workflow, switch to the `v1` branch and manage its dependencies there.

Show the new CLI:

```bash
pixi run virusdetect -- --help
```

Inspect database resolution:

```bash
pixi run virusdetect -- db path
pixi run virusdetect -- db verify
```

Check external tools on `PATH`:

```bash
pixi run virusdetect -- tools check
pixi run virusdetect -- tools install-hint
```

On `main`, `virusdetect tools check` now reports only tools found from the active environment, so missing dependencies fail fast instead of silently falling back to legacy bundled binaries.

If you need the historical Perl dependency set during migration, use the `v1` branch instead of `main`.

Bioconda packaging scaffolding now lives in [BIOCONDA.md](/Users/kentnf/projects/cornell/VirusDetect/BIOCONDA.md) and [recipes/virusdetect/meta.yaml](/Users/kentnf/projects/cornell/VirusDetect/recipes/virusdetect/meta.yaml). The public package is not published there yet, but the repository now carries the recipe and release-asset layout needed for submission.

### Alternative: local editable install

```bash
python -m pip install -e .
virusdetect --help
```

## Database Strategy

VirusDetect v2 is moving away from committing large reference databases directly into the code repository.

Planned behavior:

- database archives will be distributed separately from source code
- users will install them with `virusdetect db download`
- maintainers can generate release-ready bundles with `virusdetect db bundle`
- runtime lookup will use:
  - `VIRUSDETECT_DB_DIR` if set
  - `./database` when it already contains a valid v2 database
  - `$CONDA_PREFIX/share/virusdetect/database`
  - `~/.local/share/virusdetect/database`

The historical `databases/` directory has been removed from `main`. Use an installed database directory or a separately staged legacy database path when you need verification, bundle generation, or migration checks.

## Current CLI

```bash
virusdetect --version
virusdetect db path
virusdetect db verify
virusdetect db download
virusdetect db download --release-tag v2.0.0a0 --db-version 2026.04
virusdetect db download --url <archive-url>
virusdetect db bundle --path <db_dir> --output dist/virusdetect-db-2.0.0a1.tar.gz
virusdetect tools check
virusdetect tools install-hint
virusdetect run <reads.fa> --check-only
virusdetect run <reads.fa> --prepare-only --read-length 21-23
virusdetect run <reads.fa> --stop-after align_to_virus_reference
virusdetect run <reads.fa> --stop-after host_subtraction --host-reference host.fa
virusdetect run <reads.fa> --stop-after de_novo_assembly
virusdetect run <reads.fa> --stop-after de_novo_assembly --assembler velvet
virusdetect run <reads.fa> --stop-after generate_aligned_contigs
virusdetect run <reads.fa> --stop-after combine_contigs
virusdetect run <reads.fa> --stop-after identify_virus
virusdetect run test_data --db-path <db_dir> --stop-after identify_virus
```

Maintainer example for packaging an external historical database directory as a v2 bundle:

```bash
virusdetect db bundle --path <db_dir> --db-version 2026.04
```

This writes a tarball under `dist/` with a generated `manifest.json` and a sibling `.sha256` file.

By default, `virusdetect db download` resolves to the current VirusDetect GitHub release assets and uses the matching `.sha256` file automatically. You can override the repository, tag, or asset names with CLI flags when staging a different database bundle.

The default v2 assembler is `spades`. For v1-style comparison work, you can switch to `velvet` with `--assembler velvet`.

Recommended v1-style comparison command:

```bash
virusdetect run <reads.fa> --assembler velvet --rm-dup --stop-after identify_virus
```

## Python Identify Outputs

When you stop the Python backend after `identify_virus`, the current v2 code writes a result directory like:

```text
result_<sample>/
  contig_sequences.fa
  contig_sequences.blastn.fa
  contig_sequences.blastx.fa
  contig_sequences.undetermined.fa
  blastn.references.fa
  blastx.references.fa
  blastn.html
  blastx.html
  undetermined_blast.html
  undetermined.html
  blastn.sam
  blastx.sam
  blastn.xls
  blastx.xls
```

The default Python run also removes the matching `<sample>_temp/` working directory after a successful identify stage unless you pass `--keep-temp`.

The current HTML output provides:

- `blastn.html` and `blastx.html` summarize the final reference groups that remain after coverage/depth and redundancy filtering.
- Summary pages include compact per-reference coverage strips and stable row anchors.
- `undetermined.html` separates undetermined contigs into:
  - contigs with virus-database matches that were not kept in the final sets
  - contigs with no virus-database matches
- For small-RNA runs, undetermined reports include read-length columns (`18` to `33`), `21-22 nt` enrichment highlighting, candidate-focus sections, and stable anchors.
- `undetermined_blast.html` is the match-only view of the same undetermined set.

What is still missing from the Perl workflow:

- graphical alignment views from the legacy `Bio::Graphics` reports
- full legacy-format siRNA size-distribution reporting
- complete feature parity for every legacy reporting artifact

Example smoke test for the current Python path:

```bash
export VIRUSDETECT_SMOKE_DB_PATH=/path/to/database
pixi run smoke-identify
```

For local maintainer tasks, `pixi run package-db` and `pixi run smoke-identify` now require explicit `VIRUSDETECT_PACKAGE_DB_PATH` and `VIRUSDETECT_SMOKE_DB_PATH` environment variables so they no longer depend on the repository `databases/` directory.

See [CHANGELOG.md](/Users/kentnf/projects/cornell/VirusDetect/CHANGELOG.md) for the current alpha release summary.
See [RELEASE.md](/Users/kentnf/projects/cornell/VirusDetect/RELEASE.md) for the release checklist and tagging workflow.

Representative outputs from that smoke test include:

- `vd_identify_release/result_test_data/blastn.html`
- `vd_identify_release/result_test_data/blastx.html`
- `vd_identify_release/result_test_data/undetermined.html`
- `vd_identify_release/result_test_data/undetermined_blast.html`

## Repository Layout

```text
src/virusdetect/   Python package for v2
tests/             Unit tests for the Python rewrite line
```

## Legacy Usage

Legacy compatibility notes and historical Perl entrypoints now live on the `v1` branch.
