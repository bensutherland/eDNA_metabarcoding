#!/bin/bash
# Remove artefactual reads using obigrep based on sizes 

# Global variables
SAMPLE_FOLDER="04b_annotated_samples"

# Run on all uniq.fa files 
ls -1 $SAMPLE_FOLDER/*merged_data_assi_230_uniq_c2_55-300.fa | \
    perl -pe 's/\.fa//' | \
    sort -u | \
    while read i
    do
        echo $i
        #export
        obitab --output-seq $i".fa" > $i".txt"
    done

