#!/bin/bash

module load kraken2/2.1.6

DB="/scratch/$USER/kraken2_16gb_db"
FASTQ_DIR="/scratch/$USER/shotgun_metagenomics_A3/fastq"
OUT_DIR="/scratch/$USER/kraken_results_16gb"

mkdir -p $OUT_DIR

for sample in SRR8146935 SRR8146936 SRR8146938 SRR8146951 SRR8146952 SRR8146954
do
	echo "Processing $sample..."

	if [ -f $OUT_DIR/${sample}.report ]; then
		echo "$sample already done, skipping..."
		continue
	fi

	rm -f "$OUT_DIR/${sample}.kraken" "$OUT_DIR/${sample}.report"

	kraken2 \
        	--db $DB \
        	--paired \
        	$FASTQ_DIR/${sample}_1.fastq.gz \
        	$FASTQ_DIR/${sample}_2.fastq.gz \
        	--confidence 0.10 \
		--threads 2 \
        	--report $OUT_DIR/${sample}.report \
        	--output $OUT_DIR/${sample}.kraken

done

echo "All samples processed!"
