# Bioconda PR Draft

Title:

```text
Add virusdetect 2.0.0a0
```

Body:

```text
This PR adds a new Bioconda recipe for VirusDetect v2.

Summary:
- package: virusdetect
- version: 2.0.0a0
- source: https://github.com/kentnf/VirusDetect/releases/download/v2.0.0a0/virusdetect-2.0.0a0.tar.gz
- sha256: d741cbb5150adae0e61962f5d26c0fbecae11b10a0744c24a9558e22c3b9e4c6

Runtime dependencies:
- bwa
- samtools
- blast
- hisat2
- spades
- velvet

Notes:
- This is the Python rewrite line of VirusDetect.
- The packaged database is distributed separately from the code package and installed with `virusdetect db download`.
- `velvet` is included because the CLI exposes `--assembler velvet` as a supported alternate assembler backend.

Local validation performed:
- conda render recipes/virusdetect/meta.yaml --override-channels -c conda-forge -c bioconda
- conda build recipes/virusdetect/meta.yaml --override-channels -c conda-forge -c bioconda --no-anaconda-upload
```
