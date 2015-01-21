#!/bin/sh

# unpack plant virus database and create index
gunzip databases/vrl_plant.gz
./bin/samtools faidx databases/vrl_plant
./bin/bwa index databases/vrl_plant
./bin/formatdb -i databases/vrl_plant -p F

gunzip test_data.gz
