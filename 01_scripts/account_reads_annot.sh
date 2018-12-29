#!/bin/bash
# Count number of reads per amplicon type for 18S data in 04_samples 
# For dino
for i in 04_samples/*assi_dino_sannot.fq ; do grep -cE '^\+$' $i ; done > 07_results/assigned_dino_fastq_read_count.txt

# For diat
for i in 04_samples/*assi_diat_sannot.fq ; do grep -cE '^\+$' $i ; done > 07_results/assigned_diat_fastq_read_count.txt

