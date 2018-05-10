#!/bin/bash
# Remove the primer sequences from the forward and reverse reads 
# have not yet applied quality filter, as this is before read merging

# Global variables
RAW_FOLDER="02_raw_data"

# Define adapters (todo: use separate file and queue as file)
# for rbcL
#ADAPT1="CCRTTYATGCGTTGGAGAGA"
#ADAPT2="AARCAACCTTGTGTAAGTCT"

# for 16S
#ADAPT1="CCTACGGGNGGCWGCAG"
#ADAPT2="GACTACHVGGGTATCTAATCC"

# for 18S
#ADAPT1="AACCTGGTTGATCCTGCCAGT"
#ADAPT2="GCTATTGGAGCTGGAATTAC"

# for LSU
ADAPT1="AMAAGTACCRYGAGGGAAAG"
ADAPT2="SCWCTAATCATTCGCTTTACC"

# Merging multiple files consecutively
ls -1 $RAW_FOLDER/*.fastq.gz | \
    # Remove path
    awk -F"/" '{ print $2 }' | \
    
    # Remove end of name   
    perl -pe 's/R[12]\_001\.fastq\.gz//' | \
    sort -u | \
    
    # Run cutadapt on each input file
    # Note that primers are being removed from the 5' end of forward and reverse reads
    while read i
    do
        echo $i
        cutadapt -g $ADAPT1 -G $ADAPT2 -o $RAW_FOLDER/$i"R1_001_noprime.fastq" -p $RAW_FOLDER/$i"R2_001_noprime.fastq" $RAW_FOLDER/$i"R1_001.fastq.gz" $RAW_FOLDER/$i"R2_001.fastq.gz" 2>&1 >> 10_log_files/"cutadapt.log" 
    done

