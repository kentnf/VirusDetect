# VirusDetect v2

Python rewrite of the VirusDetect plant virus discovery pipeline.

Current package version: `2.0.0a0`

Legacy Perl releases remain on the `v1` branch and historical `v1.2` to `v1.8.1` tags. The `main` branch is the Python v2 line.

## Install

### Recommended: pixi

```bash
pixi install
pixi run virusdetect -- --version
pixi run virusdetect -- tools check
```

The locked pixi environment currently installs:

- `bwa`
- `samtools`
- `blast`
- `hisat2`
- `spades`
- `velvet`

### Alternative: editable Python install

```bash
python -m pip install -e .
virusdetect --version
```

### Bioconda status

A Bioconda recipe scaffold is kept in [recipes/virusdetect/meta.yaml](/Users/kentnf/projects/cornell/VirusDetect/recipes/virusdetect/meta.yaml). It is not published yet. Submission notes are in [BIOCONDA.md](/Users/kentnf/projects/cornell/VirusDetect/BIOCONDA.md).

## Database

VirusDetect v2 uses a separately distributed database bundle instead of storing the reference database in this repository.

Check the default install target:

```bash
virusdetect db path --target
```

Download the current packaged database:

```bash
virusdetect db download
virusdetect db verify
```

The default download source is the matching GitHub release asset for the current package version:

- release tag: `v2.0.0a0`
- database bundle: `virusdetect-db-2026.04.tar.gz`

Database lookup order at runtime:

1. `VIRUSDETECT_DB_DIR`
2. `./database`
3. `$CONDA_PREFIX/share/virusdetect/database`
4. `~/.local/share/virusdetect/database`

Maintainers can build a release bundle from an existing database directory with:

```bash
virusdetect db bundle --path <db_dir> --db-version 2026.04
```

## Quick Start

Minimal end-to-end example with the repository test dataset:

```bash
pixi run virusdetect -- db download --path ./database
pixi run virusdetect -- run test_data --db-path ./database --stop-after identify_virus -o vd_identify_release
```

The default v2 assembler is `spades`. For v1-style comparison work, use:

```bash
virusdetect run <reads.fa> --assembler velvet --rm-dup --stop-after identify_virus
```

## Result Files

Stopping after `identify_virus` produces a Perl-style release directory:

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

After a successful identify step, the matching `<sample>_temp/` directory is removed unless `--keep-temp` is used.

## Current Scope

Implemented on `main`:

- Python CLI with `run`, `db`, and `tools`
- sample preparation and validation
- virus-reference alignment
- host subtraction
- de novo assembly with `spades` or `velvet`
- aligned-contig generation
- combined-contig redundancy reduction
- Python-native `blastn` and `blastx` virus identification
- Perl-style release output set for the identify stage

Still not ported from the legacy Perl workflow:

- `Bio::Graphics` report parity
- full historical siRNA table and layout parity
- exact parity for every legacy HTML formatting detail

## Repository Layout

```text
src/virusdetect/   Python package
tests/             Unit tests
tools/             Restored standalone legacy utilities
recipes/           Bioconda recipe scaffold
```

See [CHANGELOG.md](/Users/kentnf/projects/cornell/VirusDetect/CHANGELOG.md) for release notes and [RELEASE.md](/Users/kentnf/projects/cornell/VirusDetect/RELEASE.md) for the release checklist.
