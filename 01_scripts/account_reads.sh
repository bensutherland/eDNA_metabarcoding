#!/bin/bash
# Count number of reads per sample in 02_raw_data 
for i in 02_raw_data/*R1_001.fastq ; do grep -cE '^\+$' $i ; done > 07_results/input_fastq_read_count.txt

