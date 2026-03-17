#!/bin/bash

module load fastqc

FASTQ_DIR=/scratch/$USER/shotgun_metagenomics_A3/fastq

OUTDIR=results/fastqc

mkdir -p $OUTDIR


for file in $FASTQ_DIR/*.fastq.gz
do
fastqc $file -o $OUTDIR
done
