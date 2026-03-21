#!/bin/bash

module load fastp

FASTQ_DIR=/scratch/$USER/shotgun_metagenomics_A3/fastq
OUTDIR=/scratch/$USER/shotgun_metagenomics_A3/trimmed_fastq

mkdir -p $OUTDIR

for sample in SRR8146935 SRR8146936 SRR8146938 SRR8146951 SRR8146952 SRR8146954
do
    fastp -i ${FASTQ_DIR}/${sample}_1.fastq.gz \
          -I ${FASTQ_DIR}/${sample}_2.fastq.gz \
          -o ${OUTDIR}/${sample}_trimmed_1.fastq.gz \
          -O ${OUTDIR}/${sample}_trimmed_2.fastq.gz \
          --cut_tail \
          --cut_mean_quality 20 \
          --length_required 50 \
          --html ${OUTDIR}/${sample}_fastp.html \
          --json ${OUTDIR}/${sample}_fastp.json
done
