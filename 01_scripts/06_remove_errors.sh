#!/bin/bash
# Remove artefactual reads using obigrep based on sizes 

# Global variables
SAMPLE_FOLDER="04_samples"

# Run on all uniq.fa files 
ls -1 $SAMPLE_FOLDER/*uniq_c*.fa | \
    perl -pe 's/\.fa//' | \
    sort -u | \
    while read i
    do
        echo $i
        # Label with obiclean
        obiclean -s merged_sample -r 0.05 $i".fa" > $i"clean.fa"
        # Filter
        obigrep -a 'obiclean_status:s|h' $i"clean.fa" > $i"clean_HS.fa"
    done

