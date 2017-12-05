#!/bin/bash
# Remove artefactual reads using obigrep based on sizes 

# Global variables
SAMPLE_FOLDER="04_samples"

# User defined variables
LMIN=55 # lower size limit
LMAX=75 # upper size limit
MIN_READS="10" # minimum count per sample 

# Run on all uniq.fa files 
ls -1 $SAMPLE_FOLDER/*uniq.fa | \
    perl -pe 's/\.fa//' | \
    sort -u | \
    while read i
    do
        echo $i
        obigrep --lmin $LMIN --lmax $LMAX -p 'count>='"$MIN_READS" $i".fa" > "$i"_c"$MIN_READS"_"$LMIN"-"$LMAX".fa
    done


# Move files to merged folder
#mv $MERGED_FOLDER/*ali.fq $OUTPUT_FOLDER
