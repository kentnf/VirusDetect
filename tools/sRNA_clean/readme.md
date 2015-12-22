
clean sRNA adapter 
==================

The usage of sRNA_clean.pl is quite simple, just provide adapter sequences and fastq file
```
perl sRNA_clean.pl -s CAGATCGGAAGAGCACA input.fastq
```

We collect 6 sRNA adapters from various sRNA read files (below), and only the first 17nt was taken for adapter removing.

```
CTGTAGGCACCATCAAT
CAGATCGGAAGAGCACA
TCGTATGCCGTCTTCTG
TGGAATTCTCGGGTGCC
ATCTCGTATGCCGTCTT
GTACCTCGTATGCCGTC
```

The output files include cleaned reads after remove adapter, and report file including statistics information about adapter removing. Below is the format of report file:

sample | total | unmatch | null | match | baseN | short | clean
--- | --- | --- | --- | --- | --- | --- | ---
test.fq | 1000 | 132 | 0 | 868 | 0 | 6 | 862

- **sample**: sample name
- **total**: total number of reads before adapter removing
- **unmatch**: number of reads do not contain provide adapter
- **null**: number of reads do not contain sRNA 
- **match**: number of reads contain both sRNA and adapter sequence
- **baseN**: number of match reads contain undetermined base in sRNA
- **short**: number of match reads contain short sRNA after adapter removing
- **clean**: number of match reads contain properly sRNA

