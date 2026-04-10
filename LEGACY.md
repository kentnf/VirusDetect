# Legacy Compatibility

This document tracks the historical Perl workflow after the Python rewrite became the only supported CLI path on `main`.

For long-term maintenance of the historical implementation, prefer the `v1` branch and the `v1.x` tags.

The repository default branch is now `main`. The old `master` branch is only a frozen historical alias and should not receive new work.

For the remaining code-level removal plan, see [LEGACY_REMOVAL.md](/Users/kentnf/projects/cornell/VirusDetect/LEGACY_REMOVAL.md).

The historical Perl assets have now been removed from `main`. For long-term maintenance or reproducible reruns of the Perl implementation, use the `v1` branch and the `v1.x` tags.

## Legacy Entry Points

The historical Perl commands were:

```bash
perl virus_detect.pl --reference vrl_plant test_data
perl bin/virus_identify.pl --reference databases/vrl_plant input_contig.fasta
```

Those command entrypoints have now been removed from `main`. Run them from the `v1` branch if you still need the historical Perl workflow.

## Legacy Prerequisites

The historical Perl toolchain is no longer part of the supported `main` environment.

Install the legacy extras explicitly when needed:

```bash
pixi add perl perl-bioperl
virusdetect tools install-hint --legacy
```

`Bio::Graphics` may still need extra platform-specific setup beyond `perl-bioperl`.

## Transitional Database Usage

The repository `databases/` directory has been removed from `main`.

When you intentionally want to inspect or package an external legacy-style database directory on `main`, use explicit database commands:

```bash
virusdetect run test_data --db-path <db_dir> --stop-after identify_virus
virusdetect db verify --path <db_dir>
virusdetect db bundle --path <db_dir> --db-version 2026.04
```

For supported legacy execution, use the `v1` branch. The Perl files themselves are no longer present on `main`.
