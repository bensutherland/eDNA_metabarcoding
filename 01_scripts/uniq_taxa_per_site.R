# Uses output from Rscript 'read_counts_to_annotation.R' to consider the unique taxa observed at each site

#rm(list=ls())

par(mfrow=c(2,3), mar= c(4,3,3,1) + 0.2, mgp = c(2,0.75,0))
options(scipen = 9999999)

# Choose dataset
datatypes <- c("C3_val", "C3_COI", "C3_16s",  "SOG_val", "SOG_16s")

# Collect information of uniq taxa and total reads per site
# Also plots uniq taxa by read count

# Set Nulls
taxa.read.counts.list <- list(); uniq.taxa <- NULL; reads.per.site <- NULL

for(d in datatypes){
  # Set working directory
  working.dir <- paste("~/Documents/03_eDNA/eDNA_metabarcoding_", d, sep = "")
  setwd(working.dir)
  
  # Import filename for the 'x count' filtered count file
  filename <- paste("05_annotated/", d, "_count_by_taxa_filt_at_10.csv", sep = "")
  counts <- read.delim2(filename, sep = ",", row.names = 1)
  
  # Record number of taxa per site
  uniq.taxa <- colSums(counts != 0) # if the value is 0 it is FALSE, if not 0, is TRUE, then the num of true is summed
  
  # Record number of reads per site
  reads.per.site <- colSums(counts)
  
  taxa.read.counts.list[[d]] <- rbind(uniq.taxa, reads.per.site)
  
  # plot but without mock sample
  plot(uniq.taxa[uniq.taxa != "sample.Mock"] 
       ~ reads.per.site[reads.per.site != "sample.Mock"]
       , main = d
       , xlab = "Number Reads", ylab = "Number Taxa", las = 1)
 
}


### Plot unique taxa per location for each amplicon (in progress)

# make filename this way: filename <- paste("05_annotated/", datatype, "_unassigned_unknown_counts.csv", sep = "")
# pdf(file = filename, width = 10, height = 8)

## plot taxa C3
c3.amplicons <- c("C3_16s", "C3_val", "C3_COI")
par(mfrow=c(3,1), mar=c(5,5,3,1))

for(a in c3.amplicons){
  barplot(taxa.read.counts.list[[a]][1,], las = 2, main = a
          , ylab = "Number unique taxa")
}


## plot taxa SOG
# problem with taxa naming
colnames(taxa.read.counts.list[["SOG_val"]]) <- gsub(colnames(taxa.read.counts.list[["SOG_val"]]), pattern = "sample.", replacement = "")
data <- taxa.read.counts.list[["SOG_val"]]
data <- t(data)
data <- data[order(rownames(data)),]

data <- t(data)

taxa.read.counts.list[["SOG_val"]] <- data

# For SOG_16s
taxa.read.counts.list[["SOG_16s"]]
colnames(taxa.read.counts.list[["SOG_16s"]]) <- gsub(colnames(taxa.read.counts.list[["SOG_16s"]]), pattern = "sample.S_", replacement = "")
data <- taxa.read.counts.list[["SOG_16s"]]
data <- t(data)
data <- data[order(rownames(data)),]

data <- t(data)

taxa.read.counts.list[["SOG_16s"]] <- data


# Choose amplicons
sog.amplicons <- c("SOG_val", "SOG_16s")

par(mfrow=c(2,1), mar=c(5,5,3,1))

for(a in sog.amplicons){
  barplot(taxa.read.counts.list[[a]][1,], las = 2, main = a)
  
}

# ## plot taxa all
# datatypes
# par(mfrow=c(5,1), mar=c(5,5,3,1))
# 
# # intro
# 
# for(a in datatypes){
#   barplot(taxa.read.counts.list[[a]][1,], las = 2, main = a)
#   }
# 

## write out
