# Calculate summary statistics for eDNA analyses
# Depends on previously running 01_scripts/check_merging.sh, works on final output from this script


# Clear space
# rm(list=ls())

## Install and load packages
# NA

# Set working directory
setwd("~/Documents/04_HAB/eDNA_metabarcoding_HAB_LSU_ngsfilter")

# Read in data
my.data <- read.csv2(file = "07_results/proportion_aligned_full.csv", header = F, sep = ",")
my.data <- my.data[,-1] # drop first col
head(my.data)

# Give column names
colnames(my.data) <- c("num.raw.reads", "filename", "num.merged.reads", "percent.merged")
head(my.data)
my.data$percent.merged <- as.numeric(as.character(my.data$percent.merged))

#### Summarize Raw Data ####
sum(my.data$num.raw.reads)
summary(my.data$num.raw)
sd(my.data$num.raw.reads)

#### Summarize Merged Data ####
sum(my.data$num.merged.reads)
summary(my.data$percent.merged)
sd(my.data$percent.merged)


#### Percent Merged Reads ####
sum(my.data$num.merged.reads) / sum(my.data$num.raw.reads) * 100
