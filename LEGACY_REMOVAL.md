# Legacy Removal Checklist

This file tracks the migration status of legacy Perl assets on `main`.

## Already Isolated on `main`

- `virusdetect run` on `main` is now Python-only; the CLI rejects `--backend legacy` and points users to `v1`.
- Python runs no longer import or dispatch into the legacy backend module.
- `src/virusdetect/legacy.py` and its dedicated tests have been removed from `main`.
- the top-level Perl entrypoint and Perl helper modules under `bin/` have been removed from `main`.
- the historical helper scripts under `tools/` are kept on `main` as standalone upstream/downstream utilities and are not part of the supported runtime path.
- the bundled legacy binaries under `bin/` have been removed from `main`.
- `virusdetect tools install-hint` now prints Python-only dependency guidance on `main`.
- Perl-compatible `.blastn.xls` / `.blastx.xls` compatibility exports have been restored in the default Python identify outputs on `main`.
- Perl-compatible `.blastn.sam` / `.blastx.sam` compatibility exports have been restored in the default Python identify outputs on `main`.
- Python-mode analysis plans no longer carry the `missing_legacy_perl_modules` field.
- legacy execution is no longer part of the supported `virusdetect` CLI surface on `main`.
- maintainer `pixi` tasks now require explicit database-path environment variables instead of hardcoding the repository `databases/` directory.

## Intentional Compatibility Kept on `main`

- historical standalone utilities under `tools/` are kept for reference and separate maintenance; the Python CLI does not invoke them as part of the supported `main` runtime.
- old database-layout detection in `src/virusdetect/db.py` is still kept so `main` can verify, bundle, and inspect historical database directories during migration.
- this remaining database compatibility layer is Python code for transition support, not a retained Perl runtime.

## Remaining Migration Work

### Branch and Release Management

- Keep `main` as the only active v2 development branch.
- Keep `v1` as the only supported maintenance branch for the historical Perl implementation.
- Keep `master` frozen as a temporary historical alias until downstream users have clearly moved to `main` or `v1`.
- Apply branch protection on GitHub for `main` and `v1`; avoid any new commits on `master`.
- Continue publishing `v1.x` tags only for critical legacy fixes and reserve all new feature work for `main`.

### Python Feature Parity

- Port the historical `Bio::Graphics` alignment views into a Python-supported reporting path.
- Port siRNA size-distribution tables and the related highlighting logic.
- Review historical HTML artifacts and decide which ones are still worth carrying forward in v2.
- Keep checking Python output against representative v1 runs so coverage/depth filtering stays behaviorally compatible where intended.
- Track identify/report-stage output differences in [REPORT_PARITY.md](/Users/kentnf/projects/cornell/VirusDetect/REPORT_PARITY.md).

### Database and Distribution

- Keep release assets and external database bundles documented clearly enough that `v1` users can still recover the historical workflow.
- Finalize the public release flow for v2 source packages, wheels, and database bundles.
- Complete the Bioconda publication path once the release artifact layout is stable.
- Ensure `virusdetect db download` points to reproducible release assets for both code and database updates.

### Verification and Migration Safety

- Expand regression tests around identify/report outputs as more legacy behavior is ported.
- Keep smoke-test inputs and expected report artifacts stable enough to catch report regressions between releases.
- Document any intentional v1/v2 output differences so users can distinguish missing parity from deliberate redesign.
