# VirusDetect v2

Current package version: `2.0.0a0`

This branch is the Python rewrite line for VirusDetect.

Legacy Perl releases remain available on the `v1` branch and the historical tags (`v1.2` through `v1.8.1`). The old entrypoints such as `virus_detect.pl` and `bin/virus_identify.pl` are still present in this repository during the migration, but new development should target the Python package under `src/virusdetect`.

## Status

The v2 branch currently provides:

- a Python package skeleton
- a `virusdetect` command-line entrypoint
- database path, verification, and archive download commands
- external tool discovery checks
- a `run` command that can validate inputs, classify samples, and prepare sample working files
- Python implementations of sample preparation, virus-reference alignment, host subtraction, de novo assembly, aligned contig generation, combined-contig deduplication, and virus identification
- Python-native `blastn`/`blastx` identification outputs, including raw tables, top-hit tables, contig classification FASTA files, per-sample summary TSV/JSON, HTML summary pages, per-reference detail pages, and undetermined-contig reports
- a transitional `legacy` backend that runs the existing Perl pipeline from the Python CLI

The legacy backend remains available during the migration, but the Python backend can now run end-to-end through `identify_virus`.

This is still an alpha release line: the Python pipeline is functional through identification and reporting, but report parity and packaged database delivery are not finished yet.

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
  - summary TSV/JSON
  - reference FASTA export
  - summary HTML with compact coverage maps
  - per-reference detail HTML with coordinate coverage maps
  - undetermined contig HTML reports

Not yet ported from the legacy Perl workflow:

- `Bio::Graphics` alignment images and detailed graphical layouts
- siRNA size-distribution tables and highlighting logic
- full parity for every historical HTML artifact and formatting detail
- final packaging of downloadable v2 databases separate from the repository

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
- transitional legacy Perl packages: `perl`, `perl-bioperl`

`Bio::Graphics` is still only a legacy Perl requirement and may need extra platform-specific setup outside the default `pixi` environment, especially on `osx-arm64`.

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

During the migration, `virusdetect tools check` will also fall back to bundled legacy tools in `bin/` when needed for transitional development.

When you use `--backend legacy`, the original Perl dependencies still apply. In particular, BioPerl modules such as `Bio::SeqIO` and `Bio::Graphics` must be installed.

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
  - `./databases` when it already contains a valid legacy database
  - `$CONDA_PREFIX/share/virusdetect/database`
  - `~/.local/share/virusdetect/database`

At the moment the bundled legacy `databases/` directory is still recognized by `virusdetect db path` and `virusdetect db verify` so the branch remains usable during the transition.

## Current CLI

```bash
virusdetect --version
virusdetect db path
virusdetect db verify
virusdetect db download
virusdetect db download --release-tag v2.0.0a0 --db-version 2026.04
virusdetect db download --url <archive-url>
virusdetect db bundle --path databases --output dist/virusdetect-db-2.0.0a1.tar.gz
virusdetect tools check
virusdetect tools install-hint
virusdetect run <reads.fa> --check-only
virusdetect run <reads.fa> --prepare-only --read-length 21-23
virusdetect run <reads.fa> --backend python --stop-after align_to_virus_reference
virusdetect run <reads.fa> --backend python --stop-after host_subtraction --host-reference host.fa
virusdetect run <reads.fa> --backend python --stop-after de_novo_assembly
virusdetect run <reads.fa> --backend python --stop-after generate_aligned_contigs
virusdetect run <reads.fa> --backend python --stop-after combine_contigs
virusdetect run <reads.fa> --backend python --stop-after identify_virus
virusdetect run <reads.fa> --backend legacy -o output_dir
```

Maintainer example for packaging the current legacy database as a v2 bundle:

```bash
virusdetect db bundle --path databases --db-version 2026.04
```

This writes a tarball under `dist/` with a generated `manifest.json` and a sibling `.sha256` file.

By default, `virusdetect db download` resolves to the current VirusDetect GitHub release assets and uses the matching `.sha256` file automatically. You can override the repository, tag, or asset names with CLI flags when staging a different database bundle.

## Python Identify Outputs

When you stop the Python backend after `identify_virus`, the current v2 code writes a result directory like:

```text
result_<sample>/
  <sample>.summary.tsv
  <sample>.summary.json
  blastn.html
  blastx.html
  undetermined.html
  undetermined_blast.html
  blastn.reference.fa
  blastx.reference.fa
  blastn_references/
  blastx_references/
```

The current HTML output provides:

- `blastn.html` and `blastx.html` summarize the final kept references after coverage/depth and redundancy filtering
- summary pages now include compact per-reference coverage strips so you can see fragmented versus near-complete support at a glance
- `blastn_references/` and `blastx_references/` contain one detail page per kept reference with coordinate coverage maps, strand-aware hit blocks, and contig tables
- `undetermined.html` separates undetermined contigs into:
  - contigs with virus-database hits that were not retained in the final sets
  - contigs with no retained virus-database hits
- `undetermined_blast.html` is the hits-only view of the same undetermined set

What is still missing from the Perl workflow:

- graphical alignment views from the legacy `Bio::Graphics` reports
- full siRNA size-distribution reporting
- complete feature parity for every legacy reporting artifact

Example smoke test for the current Python path:

```bash
pixi run virusdetect -- run test_data --backend python --stop-after identify_virus -o vd_identify_py
```

See [CHANGELOG.md](/Users/kentnf/projects/cornell/VirusDetect/CHANGELOG.md) for the current alpha release summary.
See [RELEASE.md](/Users/kentnf/projects/cornell/VirusDetect/RELEASE.md) for the release checklist and tagging workflow.

Representative outputs from that smoke test include:

- `vd_identify_py/result_test_data/blastn.html`
- `vd_identify_py/result_test_data/blastx.html`
- `vd_identify_py/result_test_data/blastn_references/AB083196.html`
- `vd_identify_py/result_test_data/blastx_references/AEV92950.1.html`
- `vd_identify_py/result_test_data/undetermined.html`
- `vd_identify_py/result_test_data/undetermined_blast.html`

## Repository Layout

```text
src/virusdetect/   Python package for v2
tests/             Unit tests for the Python CLI scaffold
bin/               Legacy bundled binaries used by v1
databases/         Legacy bundled reference data used by v1
tools/             Legacy helper scripts for preprocessing and DB generation
```

## Legacy Usage

The legacy Perl workflow is still available on this branch while the migration is in progress:

```bash
perl virus_detect.pl --reference vrl_plant test_data
perl bin/virus_identify.pl --reference databases/vrl_plant input_contig.fasta
```
