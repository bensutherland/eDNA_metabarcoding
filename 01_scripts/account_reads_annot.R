# Clear space
# rm(list=ls())

## Install and load packages
# NA

# Set working directory
setwd("~/Documents/04_HAB/eDNA_metabarcoding_HAB_18S_ngsfilter")

# Read in data
my.data <- read.csv2(file = "07_results/assigned_dino_fastq_read_count.txt", header = F, sep = ",")
my.data <- read.csv2(file = "07_results/assigned_diat_fastq_read_count.txt", header = F, sep = ",")

head(my.data)

# Summarize
sum(my.data$V1)
summary(my.data$V1)
sd(my.data$V1)
