#!/bin/bash
# Using the labels 'dino | diat' split the R1 fastq file into separate files

# Global variables
SAMPLE_FOLDER="04_samples"

# Merging multiple files consecutively
ls -1 $SAMPLE_FOLDER/*_R1_assi.fq | \
    # Remove path
    awk -F"/" '{ print $2 }' | \
    
    # Remove end of name   
    perl -pe 's/R1\_assi\.fq//' | \
    sort -u | \
    
    # Sort into dino or diat data 
    while read i
    do
        echo $i

        # Take all dino and put in one file
        grep -E -A3 'sample\=dino' $SAMPLE_FOLDER/$i"R1_assi.fq" | grep -vE '^--$' - > $SAMPLE_FOLDER/$i"R1_assi_dino.fq"

        # Take all diat and put in second file
        grep -E -A3 'sample\=diat' $SAMPLE_FOLDER/$i"R1_assi.fq" | grep -vE '^--$' - > $SAMPLE_FOLDER/$i"R1_assi_diat.fq"

    done

