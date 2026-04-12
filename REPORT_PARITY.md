# Report Parity Checklist

This document tracks identify/report-stage parity between the historical Perl `v1` line and the Python implementation on `main`.

It is intentionally limited to reporting behavior and output artifacts. Database migration support is tracked separately in [LEGACY_REMOVAL.md](/Users/kentnf/Desktop/projects/cornell/VirusDetect/LEGACY_REMOVAL.md).

## Current Coverage

The Python identify stage currently writes these primary report artifacts:

- `contig_sequences.fa`
- `contig_sequences.blastn.fa`
- `contig_sequences.blastx.fa`
- `contig_sequences.undetermined.fa`
- `blastn.references.fa`
- `blastx.references.fa`
- `blastn.html`
- `blastx.html`
- `undetermined_blast.html`
- `undetermined.html`
- `blastn.sam`
- `blastx.sam`
- `blastn.xls`
- `blastx.xls`

## Already Matched Or Replaced

- Summary pages now provide stable per-reference row anchors and compact coverage maps without depending on extra per-reference output directories.
- Undetermined-contig pages now provide grouped headers, browser-native tooltips, stable anchors, candidate-focus sections, and read-length distribution columns for small-RNA mode.
- The old Perl `Subject` wording has been replaced with `Reference` wording in current user-facing reports.
- The old `hit` / `no-hit` wording on current pages has been tightened to `match` / `no-match` where those pages describe current Python outputs.

## Remaining Gaps

- Exact `Bio::Graphics` look-and-feel is not reproduced; the Python SVG panels replace the old graphics but do not match them pixel-for-pixel.
- The small-RNA undetermined tables do not fully match the historical Perl layout and presentation.
- Legacy per-reference detail HTML pages are not part of the default Python release output.
- Some historical report behaviors still need a product decision: reproduce exactly, simplify, or retire.

## Verification

The current report stack should be checked with:

- `PYTHONPATH=src python3 -m unittest discover -s tests -v`
- `python3 -m compileall src`
- `pixi run smoke-identify`

## Regression Focus

When report behavior changes, prioritize checks that confirm:

- summary pages keep stable anchors without depending on extra result files
- undetermined reports still expose the expected match/no-match sections and small-RNA columns
- candidate-row highlighting and counts stay consistent
- final result directories only contain the documented release artifacts
