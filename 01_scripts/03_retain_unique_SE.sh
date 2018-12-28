#!/bin/bash
# Keep only unique reads per sample
# Currently applied using -m <TAG NAME> where TAG NAME is 'sample'


# Global variables
INTERP_FOLDER="00_archive"
SAMPLE_FOLDER="04b_annotated_samples"

# Run on all the assigned fastq files 
ls -1 $SAMPLE_FOLDER/*assi_230*.fq | \
    perl -pe 's/\.fq//' | \
    sort -u | \
    while read i
    do
        echo $i
        obiuniq -m sample $i".fq" > $i"_uniq.fa" 
    done


# Move files to merged folder
#mv $MERGED_FOLDER/*ali.fq $OUTPUT_FOLDER
