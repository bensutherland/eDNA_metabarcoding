#!/bin/bash
# Retain only aligned reads from multiple datafiles at once

# Global variables
INTERP_FOLDER="00_archive"
RAW_FOLDER="02_raw_data"
MERGED_FOLDER="03_merged"
OUTPUT_FOLDER="04_samples"

# Merging multiple files consecutively
ls -1 $RAW_FOLDER/*.fastq | \
    # Remove path
    awk -F"/" '{ print $2 }' | \
    
    # Remove end of name   
    perl -pe 's/R[12]\_001\.fastq//' | \
    sort -u | \
    
    # Run ngsfilter on each input file
    while read i
    do
        echo $i
        INTERP=$i"interp.txt"
        ngsfilter -t $INTERP_FOLDER/$INTERP -u $OUTPUT_FOLDER/$i"unidentified.fq" $MERGED_FOLDER/$i"ali.fq" > ./$OUTPUT_FOLDER/$i"ali_assi.fq" 
    done


# Move files to merged folder
#mv $MERGED_FOLDER/*ali.fq $OUTPUT_FOLDER
