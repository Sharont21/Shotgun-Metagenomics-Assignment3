#!/bin/bash

module load sra-toolkit/3.0.9

for sample in SRR8146935 SRR8146936 SRR8146938 SRR8146951 SRR8146952 SRR8146954
do
    prefetch $sample
done
