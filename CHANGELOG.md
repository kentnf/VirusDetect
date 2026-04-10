# Changelog

All notable changes to the Python v2 rewrite line will be documented in this file.

The legacy Perl release history remains on the `v1` branch and the historical `v1.x` tags.

## [Unreleased]

### Added

- Bioconda recipe scaffold under `recipes/virusdetect/meta.yaml`
- release packaging notes in `BIOCONDA.md` and a repository `LICENSE` file for downstream packaging metadata
- `virusdetect db bundle` to package a legacy or v2 database directory into a release-ready `database/` tarball
- generated database manifests and sibling `.sha256` checksum files for packaged database archives
- default GitHub release asset resolution for `virusdetect db download`, including automatic checksum URL lookup
- `pixi run package-db` maintainer task for building a database bundle from an explicit database path

### Changed

- `virusdetect run` on `main` is now Python-only; users who still need the historical Perl workflow should use the `v1` branch
- default `pixi` and Bioconda dependency sets on `main` no longer include Perl; legacy extras must be installed explicitly
- tool resolution on `main` no longer falls back to legacy bundled binaries in `bin/`; runtime tools must come from the active environment
- the deprecated `--backend legacy` override is now rejected on `main` and points users to the `v1` branch
- legacy Perl usage notes have been moved out of the main README into `LEGACY.md`
- `virusdetect run --help` on `main` now presents a Python-only interface
- `virusdetect tools install-hint` now defaults to Python-only dependency guidance; historical Perl extras move behind `--legacy`
- the Python-side legacy bridge module and its dedicated tests have been removed from `main`
- the historical Perl entrypoints and helper modules (`virus_detect.pl`, `bin/virus_identify.pl`, `bin/Util.pm`, `bin/align.pm`, `bin/Bio/Graphics.pm`) have been removed from `main`
- the repository `databases/` directory has been removed from `main`; maintainers now package external database paths explicitly
- the historical helper scripts under `tools/` and the bundled legacy binaries under `bin/` have been removed from `main`
- Python-mode analysis plans no longer include `missing_legacy_perl_modules`
- maintainer `pixi` tasks now require explicit database-path environment variables instead of hardcoding the repository `databases/` directory

## [2.0.0a0] - 2026-04-10

Initial public alpha snapshot of the Python rewrite line on `main`.

### Added

- Python package layout under `src/virusdetect`
- `virusdetect` CLI with `run`, `db`, and `tools` subcommands
- database path resolution, verification, and archive download plumbing
- staged Python pipeline execution through:
  - sample preparation
  - virus-reference alignment
  - host subtraction
  - de novo assembly
  - aligned contig generation
  - combined-contig deduplication
  - virus identification
- Python-native `blastn` and `blastx` identification outputs
- summary TSV and JSON outputs for final kept references
- reference FASTA export for final `blastn` and `blastx` sets
- HTML summary reports with compact coverage maps
- per-reference HTML detail reports with strand-aware coordinate coverage maps
- undetermined contig reports split into:
  - contigs with virus-database hits
  - contigs without retained virus-database hits
- transitional legacy backend that still invokes the original Perl workflow from the new CLI
- pixi environment definition using `conda-forge` and `bioconda`

### Changed

- version metadata is now aligned on the Python v2 line as `2.0.0a0`
- install hints now reflect that `hisat2` and `spades` are used by the current Python pipeline
- README now documents the current v2 migration status, report outputs, and known gaps

### Known Limitations

- legacy `Bio::Graphics` image reports are not yet ported
- siRNA size-distribution reporting and highlighting are not yet ported
- downloadable packaged v2 database releases are not finalized yet
- report formatting is intentionally simpler than the historical Perl output

### Validation

- unit tests: `PYTHONPATH=src python3 -m unittest discover -s tests -v`
- repeated smoke tests on `test_data` through `--backend python --stop-after identify_virus`
