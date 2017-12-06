#!/bin/bash
# Remove artefactual reads using obigrep based on sizes 

# Global variables
SAMPLE_FOLDER="04_samples"

# Run on all uniq.fa files 
ls -1 $SAMPLE_FOLDER/*clean_HS.fa | \
    perl -pe 's/\.fa//' | \
    sort -u | \
    while read i
    do
        echo $i
        #export
        obitab --output-seq $i".fa" > $i".txt"
    done

