# Release Guide

This document is the working release checklist for the Python v2 line.

Current target release: `2.0.0a0`

## Scope

`2.0.0a0` is an alpha release of the Python rewrite on `main`.

Included in this alpha:

- Python CLI and staged pipeline execution
- Python path through `identify_virus`
- summary TSV/JSON outputs
- HTML summary reports
- per-reference detail HTML reports
- undetermined contig reports
- transitional legacy backend from the new CLI

Not included in this alpha:

- legacy `Bio::Graphics` image reports
- siRNA size-distribution reporting
- public Bioconda package publication
- full report parity with the Perl workflow

## Pre-Release Checklist

1. Confirm version metadata matches the intended release:
   - [pyproject.toml](/Users/kentnf/projects/cornell/VirusDetect/pyproject.toml)
   - [src/virusdetect/__init__.py](/Users/kentnf/projects/cornell/VirusDetect/src/virusdetect/__init__.py)
2. Confirm release notes are updated:
   - [CHANGELOG.md](/Users/kentnf/projects/cornell/VirusDetect/CHANGELOG.md)
   - [README.md](/Users/kentnf/projects/cornell/VirusDetect/README.md)
3. Refresh the environment:

```bash
pixi install
```

4. Run the unit test suite:

```bash
pixi run test
```

5. Verify CLI version output:

```bash
pixi run virusdetect -- --version
```

6. Run the Python identify smoke test:

```bash
pixi run smoke-identify
```

7. Inspect the key smoke-test artifacts:

```text
vd_identify_release/result_test_data/blastn.html
vd_identify_release/result_test_data/blastx.html
vd_identify_release/result_test_data/blastn_references/AB083196.html
vd_identify_release/result_test_data/blastx_references/AEV92950.1.html
vd_identify_release/result_test_data/undetermined.html
vd_identify_release/result_test_data/undetermined_blast.html
vd_identify_release/result_test_data/test_data.summary.json
```

8. Build the source and wheel distributions:

```bash
pixi run package
```

9. Upload the Python package assets to GitHub Release:

```bash
gh release upload v2.0.0a0 dist/virusdetect-2.0.0a0.tar.gz dist/virusdetect-2.0.0a0-py3-none-any.whl --clobber
```

10. Build the database bundle when preparing a database release:

```bash
pixi run package-db
```

11. Upload the database archive and checksum to the matching GitHub release:

```bash
gh release upload v2.0.0a0 dist/virusdetect-db-2026.04.tar.gz dist/virusdetect-db-2026.04.tar.gz.sha256 --clobber
```

12. Verify the generated artifacts exist under `dist/` and on the release page.

13. Update [recipes/virusdetect/meta.yaml](/Users/kentnf/projects/cornell/VirusDetect/recipes/virusdetect/meta.yaml) and validate it with `conda render` when preparing a Bioconda submission.

14. Tag the release only after the worktree is reviewed and committed.

## Suggested Tagging

For this alpha line:

```bash
git tag -a v2.0.0a0 -m "VirusDetect Python alpha release"
```

Push the branch and tag when ready:

```bash
git push origin main
git push origin v2.0.0a0
```

## Release Notes Template

Title:

```text
VirusDetect v2.0.0a0
```

Summary:

```text
Alpha release of the Python rewrite line. This release provides a functional Python CLI and pipeline through virus identification, including HTML summary reports, per-reference detail reports, and undetermined-contig reports. Legacy Perl execution remains available through the transitional backend, while graphical Bio::Graphics reports and siRNA-size reporting are still pending.
```

Validation:

```text
- pixi run test
- pixi run smoke-identify
- pixi run package
```
