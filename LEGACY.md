# Legacy Compatibility

This document tracks the transitional Perl compatibility path that still exists on `main`.

For long-term maintenance of the historical implementation, prefer the `v1` branch and the `v1.x` tags.

## Legacy Entry Points

The old Perl commands are still present in this repository:

```bash
perl virus_detect.pl --reference vrl_plant test_data
perl bin/virus_identify.pl --reference databases/vrl_plant input_contig.fasta
```

## Legacy Prerequisites

The compatibility path is no longer part of the default `main` environment.

Install the legacy extras explicitly when needed:

```bash
pixi add perl perl-bioperl
```

`Bio::Graphics` may still need extra platform-specific setup beyond `perl-bioperl`.

## Transitional Database Usage

The repository `databases/` directory is no longer auto-discovered on `main`.

When you intentionally want to run transitional checks against it, pass it explicitly:

```bash
virusdetect run test_data --db-path databases --backend legacy
virusdetect run test_data --db-path databases --stop-after identify_virus
virusdetect db verify --path databases
```
