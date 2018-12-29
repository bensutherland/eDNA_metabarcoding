#!/bin/bash
# Retain only aligned reads from multiple datafiles at once (single-end read only)

# Global variables
INTERP_FOLDER="00_archive"
RAW_FOLDER="02_raw_data"
MERGED_FOLDER="03_merged"
OUTPUT_FOLDER="04_samples"

INTERP_FORWARD="interp_forward_18S.txt"

# Use ngsfilter on forward reads 
ls -1 $RAW_FOLDER/*R1_001.fastq | \

    # Remove path
    awk -F"/" '{ print $2 }' | \
    
    # Remove end of name   
    perl -pe 's/\_001\.fastq//' | \
    sort -u | \
    
    # Run ngsfilter on each input file
    while read i
    do
        echo $i
        ngsfilter -t $INTERP_FOLDER/$INTERP_FORWARD -u $OUTPUT_FOLDER/$i"unidentified.fq" $RAW_FOLDER/$i"_001.fastq" > ./$OUTPUT_FOLDER/$i"_assi.fq" 
    done

