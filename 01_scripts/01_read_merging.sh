#!/bin/bash
# Merging multiple files at once

# Global variables
RAW_FOLDER="02_raw_data"
OUTPUT_FOLDER="03_merged"

# Launch obitools
obitools

# Merging multiple files consecutively
ls -1 $RAW_FOLDER/*.fastq | \
    perl -pe 's/R[12]\_001\.fastq//' | \
    sort -u | \
    while read i 
    do
        echo $i
        illuminapairedend --score-min=40 \
        -r $i"R1_001.fastq" \
        $i"R2_001.fastq" > $i"merged.fq" 
    done

# Move files to merged folder
mv $RAW_FOLDER/*merged.fq $OUTPUT_FOLDER
