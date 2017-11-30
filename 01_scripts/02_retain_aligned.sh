#!/bin/bash
# Retain only aligned reads from multiple datafiles at once

# Global variables
MERGED_FOLDER="03_merged"

# Merging multiple files consecutively
ls -1 $MERGED_FOLDER/*merged.fq | \
    perl -pe 's/merged\.fq//' | \
    sort -u | \
    while read i 
    do
        echo $i
    # Use obigrep to only keep aligned sequences
        obigrep -p 'mode!="joined"' $i"merged.fq" > $i"ali.fq"
    done

# Move files to merged folder
#mv $MERGED_FOLDER/*ali.fq $OUTPUT_FOLDER
