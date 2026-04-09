# Changelog

All notable changes to the Python v2 rewrite line will be documented in this file.

The legacy Perl release history remains on the `v1` branch and the historical `v1.x` tags.

## [Unreleased]

### Added

- `virusdetect db bundle` to package a legacy or v2 database directory into a release-ready `database/` tarball
- generated database manifests and sibling `.sha256` checksum files for packaged database archives
- default GitHub release asset resolution for `virusdetect db download`, including automatic checksum URL lookup
- `pixi run package-db` maintainer task for building the default database bundle from `databases/`

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
