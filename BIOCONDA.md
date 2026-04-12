# Bioconda Submission Notes

Bioconda recipe scaffold: [recipes/virusdetect/meta.yaml](/Users/kentnf/projects/cornell/VirusDetect/recipes/virusdetect/meta.yaml)

Current release state:

- package version: `2.0.0a0`
- database version: `2026.04`
- source asset: `https://github.com/kentnf/VirusDetect/releases/download/v2.0.0a0/virusdetect-2.0.0a0.tar.gz`
- source SHA256: `d741cbb5150adae0e61962f5d26c0fbecae11b10a0744c24a9558e22c3b9e4c6`
- database asset: `https://github.com/kentnf/VirusDetect/releases/download/v2.0.0a0/virusdetect-db-2026.04.tar.gz`

## Verified Checks

The current repository state was checked with:

```bash
pixi install --locked
pixi run virusdetect -- --version
pixi run virusdetect -- tools check
pixi run virusdetect -- db path --target
pixi run virusdetect -- db download --path /tmp/virusdetect-db-check/database
pixi run virusdetect -- db verify --path /tmp/virusdetect-db-check/database
pixi run test
python -m compileall src
conda render recipes/virusdetect/meta.yaml --override-channels -c conda-forge -c bioconda
conda build recipes/virusdetect/meta.yaml --override-channels -c conda-forge -c bioconda --no-anaconda-upload
```

Observed results:

- `pixi install --locked` succeeds
- all runtime tools are available from the pixi environment
- default database download resolves to the GitHub release asset and installs a valid `manifest.json`-based database
- unit tests pass
- Python sources compile cleanly
- `conda render` passes
- `conda build` passes, including recipe test commands

## Recipe Notes

The recipe currently declares the v2 runtime stack as dependencies:

- `bwa`
- `samtools`
- `blast`
- `hisat2`
- `spades`
- `velvet`

`velvet` is included because the v2 CLI exposes `--assembler velvet` as a supported comparison backend.

The database is intentionally distributed outside the Bioconda package and installed with `virusdetect db download`.

## Release Checklist

1. Build package artifacts:

```bash
pixi run package
```

2. Build the database bundle:

```bash
export VIRUSDETECT_PACKAGE_DB_PATH=/path/to/database
pixi run package-db
```

3. Upload the source package, wheel, and database assets to the matching GitHub release:

```bash
gh release upload v2.0.0a0 \
  dist/virusdetect-2.0.0a0.tar.gz \
  dist/virusdetect-2.0.0a0-py3-none-any.whl \
  dist/virusdetect-db-2026.04.tar.gz \
  dist/virusdetect-db-2026.04.tar.gz.sha256 \
  --clobber
```

4. Refresh the recipe checksum if the source tarball changed:

```bash
shasum -a 256 dist/virusdetect-<version>.tar.gz
```

5. Validate the recipe locally when `conda-build` is available:

```bash
conda render recipes/virusdetect/meta.yaml --override-channels -c conda-forge -c bioconda
conda build recipes/virusdetect/meta.yaml --override-channels -c conda-forge -c bioconda --no-anaconda-upload
```

6. Copy both recipe files into `bioconda-recipes/recipes/virusdetect/`:

```text
recipes/virusdetect/meta.yaml
recipes/virusdetect/LICENSE
```

7. Open the Bioconda pull request.
