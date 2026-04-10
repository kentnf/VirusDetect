# Legacy Removal Checklist

This file tracks the migration status of legacy Perl assets on `main`.

## Already Isolated on `main`

- `virusdetect run` on `main` is now Python-only; the CLI rejects `--backend legacy` and points users to `v1`.
- Python runs no longer import or dispatch into the legacy backend module.
- `src/virusdetect/legacy.py` and its dedicated tests have been removed from `main`.
- the top-level Perl entrypoint and Perl helper modules under `bin/` have been removed from `main`.
- the historical helper scripts under `tools/` have been removed from `main`.
- the bundled legacy binaries under `bin/` have been removed from `main`.
- `virusdetect tools install-hint` now defaults to Python dependencies; historical Perl extras require `--legacy`.
- Python-mode analysis plans no longer carry the `missing_legacy_perl_modules` field.
- legacy execution is no longer part of the supported `virusdetect` CLI surface on `main`.
- maintainer `pixi` tasks now require explicit database-path environment variables instead of hardcoding the repository `databases/` directory.

## Remaining Migration Work

- Keep release assets and external database bundles documented clearly enough that `v1` users can still recover the historical workflow.
