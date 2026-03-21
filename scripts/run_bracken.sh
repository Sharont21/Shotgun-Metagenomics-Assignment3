#!/bin/bash

module load bracken

DB="/scratch/$USER/kraken2_16gb_db"
IN_DIR="/scratch/$USER/kraken_results_16gb"
OUT_DIR="/scratch/$USER/kraken_results_16gb"

for sample in SRR8146935 SRR8146936 SRR8146938 SRR8146951 SRR8146952 SRR8146954
do
    echo "Running Bracken on $sample..."

    bracken \
        -d $DB \
        -i $IN_DIR/${sample}.report \
        -o $OUT_DIR/${sample}.bracken \
        -r 150 \
        -l S

done

echo "Bracken complete!"
