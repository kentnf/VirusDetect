# Bioconda Packaging

This repository keeps a Bioconda-ready recipe scaffold under [recipes/virusdetect/meta.yaml](/Users/kentnf/projects/cornell/VirusDetect/recipes/virusdetect/meta.yaml).

## Current Source Asset

The recipe is pinned to the uploaded source distribution for:

- package version: `2.0.0a0`
- source asset: `https://github.com/kentnf/VirusDetect/releases/download/v2.0.0a0/virusdetect-2.0.0a0.tar.gz`
- source SHA256: `d741cbb5150adae0e61962f5d26c0fbecae11b10a0744c24a9558e22c3b9e4c6`

## Local Checklist

1. Build release assets:

```bash
pixi run package
export VIRUSDETECT_PACKAGE_DB_PATH=/path/to/database
pixi run package-db
```

2. Upload the Python package assets to the matching GitHub release:

```bash
gh release upload v2.0.0a0 \
  dist/virusdetect-2.0.0a0.tar.gz \
  dist/virusdetect-2.0.0a0-py3-none-any.whl \
  --clobber
```

3. Upload the database assets:

```bash
gh release upload v2.0.0a0 \
  dist/virusdetect-db-2026.04.tar.gz \
  dist/virusdetect-db-2026.04.tar.gz.sha256 \
  --clobber
```

4. Refresh the recipe checksum when the Python package version changes:

```bash
shasum -a 256 dist/virusdetect-<version>.tar.gz
```

5. Validate the recipe locally when `conda-build` is available:

```bash
conda render recipes/virusdetect/meta.yaml
```

6. Copy the recipe into `bioconda-recipes/recipes/virusdetect/` and submit a pull request there.

## Packaging Notes

- The recipe intentionally declares the external command-line tools as runtime dependencies so a Bioconda install pulls the current VirusDetect v2 stack in one step.
- The package remains `noarch: python`; compiled tools are provided by the runtime dependencies.
- Legacy Perl support is no longer part of the default runtime dependency set for the v2 package.
- The release database is still distributed as a separate asset and installed with `virusdetect db download`.
