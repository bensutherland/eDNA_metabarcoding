#!/bin/bash
# Denoise data based on count, amplicon size, and singleton|head|internal status 

# Global variables
SAMPLE_FOLDER="04_samples"

# User defined variables
LMIN=55 # lower size limit
LMAX=75 # upper size limit
MIN_READS="10" # minimum count per sample 

# Run on all *uniq.fa files 
ls -1 $SAMPLE_FOLDER/*uniq.fa | \
    perl -pe 's/\.fa//' | \
    sort -u | \
    while read i
    do
        echo "Treating "$i".fa" 

        # Denoise by length and depth
        echo "Denoising by length and depth"
        DENOISE_FILENAME=$i"_c"$MIN_READS"_"$LMIN"-"$LMAX".fa"
        obigrep --lmin $LMIN --lmax $LMAX \
           -p 'count>='"$MIN_READS" $i".fa" > $DENOISE_FILENAME 
        
        # Label status (singleton, head, internal) with obiclean
        echo "Removing PCR/sequencing errors"
        CLEAN_FILENAME=${DENOISE_FILENAME%.fa}"_clean.fa"
        obiclean -s merged_sample -r 0.05 $DENOISE_FILENAME > $CLEAN_FILENAME

        # Filter based on obiclean status
        HS_FILENAME=${CLEAN_FILENAME%.fa}"_HS.fa"
        obigrep -a 'obiclean_status:s|h' $CLEAN_FILENAME > $HS_FILENAME 
    done

