#!/bin/bash

# The following is to track how many reads retain through the process (for HABs)

#TODO
# 1. add the name of the directory (e.g. LSU) into the filename output
# 2. fix when gunzip there is no name of the sample
# 3. delete intermediate files
# 4. better comment below
# 5. add header to the file 

# Determine how many lines in input
for i in 02_raw_data/*R1_001.fastq.gz ; do gunzip -c $i | wc -l ; done > 07_results/raw_fastq_linecount.txt

# corrected to number of records
awk '{print $2 ", " $1 / 4 }' 07_results/raw_fastq_linecount.txt > 07_results/raw_fastq_recordcount.txt


# Determine how many reads in merged data
for i in 03_merged/*ali.fq ; do wc -l $i ; done > 07_results/ali_fq_linecount.txt
awk '{print $2 ", " $1 / 4}' 07_results/ali_fq_linecount.txt > 07_results/ali_fq_recordcount.txt


# Merge two files (note: assumes same order) 
paste -d "," 07_results/raw_fastq_recordcount.txt 07_results/ali_fq_recordcount.txt > 07_results/proportion_aligned.csv

# Calculate percentages
awk -F"," '{ print $0 "," ($4 / $2) * 100 }' 07_results/proportion_aligned.csv > 07_results/proportion_aligned_full.csv

# Then go to the script 'read_and_alignment_summary.R'
