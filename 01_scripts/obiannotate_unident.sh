#!/bin/bash
# For LSU, annotate the '_unidentified.fq' files with sample IDs in the fasta header before concatenating all together 
# The output will be $OUTPUT_FOLDER/merged_data_assi.fa

# Global variables
SAMPLE_FOLDER="04_samples"
OUTPUT_FOLDER="04b_annotated_samples"

# issue: if run twice, it will redo on top of the new files because of ambiguity in filename
# HAB-9-18S_S98_L001_R1_assi_diat.fq 
# vs
# HAB-9-18S_S98_L001_R1_assi_diat_sannot.fq

# Annotate the merged files with the sample name
ls -1 $SAMPLE_FOLDER/*_unidentified.fq | \
    perl -pe 's/\_L001\_R1\_unidentified\.fq//' | \
    
    # Remove the directory part of name 
    awk -F/ '{ print $2 }' - | \

    sort -u | \
    while read i 
    do
        echo $i

        # Run obiannotate 
        obiannotate -S sample:"$i" $SAMPLE_FOLDER/$i"_L001_R1_unidentified.fq" > $SAMPLE_FOLDER/$i"_sannot.fq"
    done

# Move files to merged folder separately for dino or diat
cat $SAMPLE_FOLDER/*_sannot.fq > $OUTPUT_FOLDER/merged_data_assi.fq

