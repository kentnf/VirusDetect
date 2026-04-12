# Changelog

All notable changes to the Python v2 rewrite line will be documented in this file.

The legacy Perl release history remains on the `v1` branch and the historical `v1.x` tags.

## [Unreleased]

### Added

- Bioconda recipe scaffolding and release-packaging notes for the Python v2 line.
- `virusdetect db bundle` plus generated manifest and checksum output for packaging external database directories.
- Default GitHub release asset resolution for `virusdetect db download`.
- `pixi run package-db` maintainer task for building a database bundle from an explicit database path.

### Changed

- `main` is now fully Python-only for supported execution; users who still need the historical Perl workflow should use the `v1` branch.
- Perl runtime guidance, bundled binaries, bridge modules, and deprecated backend hooks have been removed from `main`; standalone utilities under `tools/` are kept for separate upstream/downstream use.
- Default dependency guidance on `main` is now Python-only, and maintainer tasks now require explicit database-path environment variables.
- The repository no longer carries bundled `databases/` contents; maintainers package external database directories explicitly.
- Perl-style compatibility artifacts `blastn.sam`, `blastx.sam`, `blastn.xls`, and `blastx.xls` have been restored to the default identify output set for downstream continuity.
- Identify outputs are now restricted to the documented Perl-style release set; internal tables, raw BLAST dumps, depth tables, staging FASTA files, and per-reference detail directories are removed from the final result directory.
- Internal contig-depth intermediates (`.bwa.sai`, `.sam`, `.bam`, `.sorted.bam`, `.pileup`) are now cleaned from default identify outputs unless `--debug` is enabled.
- Reference FASTA exports now use Perl-style plural names: `blastn.references.fa` and `blastx.references.fa`.
- Per-sample `<sample>_temp/` working directories are now removed automatically after a successful identify-stage run unless `--keep-temp` is enabled.
- Summary and undetermined reports now provide stronger report parity through stable anchors, grouped headers, compact coverage maps, and clearer current-Python wording without depending on extra per-reference output directories.
- Undetermined-contig reports now use clearer `match` / `no-match` wording, grouped headers with tooltips, candidate-focus sections, and stable anchors.

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
- summary TSV and JSON outputs for final reference groups that remain after filtering
- reference FASTA export for final `blastn` and `blastx` sets
- HTML summary reports with compact coverage maps
- undetermined contig reports split into:
  - contigs with virus-database matches
  - contigs without virus-database matches
- database-layout compatibility for historical VirusDetect database directories
- pixi environment definition using `conda-forge` and `bioconda`

### Changed

- version metadata is now aligned on the Python v2 line as `2.0.0a0`
- install hints now reflect that `hisat2` and `spades` are used by the current Python pipeline
- README now documents the current v2 migration status, report outputs, and known gaps

### Known Limitations

- exact legacy `Bio::Graphics` parity is not finished
- the undetermined sRNA report now includes length-distribution columns and candidate highlighting, but full legacy-format parity is not finished
- downloadable packaged v2 database releases are not finalized yet
- report formatting is intentionally simpler than the historical Perl output

### Validation

- unit tests: `PYTHONPATH=src python3 -m unittest discover -s tests -v`
- repeated smoke tests on `test_data` through `--stop-after identify_virus`
