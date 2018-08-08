# Connect read counts to the annotations and plot
# Input: output of obitab (read counts per sample) and MEGAN (taxonomy ID per amplicon)

#rm(list=ls())

# Choose dataset
#datatype <- "C3_16s"
#datatype <- "C3_COI"
#datatype <- "SOG_16s"
#datatype <- "C3_val"
#datatype <- "SOG_val"
datatype <- "SOG_COI"

#### 0. Setup (no changes reqd) ####
# Install Packages
#install.packages("RColorBrewer")
library("RColorBrewer")

# Set working directory depending on the dataset
working.dir <- paste("/hdd/03_eDNA/eDNA_metabarcoding_", datatype, sep = "")
setwd(working.dir)
#setwd("~/Documents/03_eDNA/eDNA_metabarcoding_C3_16s_one_ended") # only req one end to assign

## Create a filenames list that contains file names for each dataset (1st obitab output; 2nd MEGAN output)
filenames.list <- list()
filenames.list[["C3_val"]] <- setNames(object = c("NGSLib6_S1_L001_ali_assi_uniq_c10_55-75_clean_HS.txt", "NGSLib6_S1_L001_ali_assi_uniq_c10_55-75clean_HS_annot-ex_sp.txt")
                                       , nm = c("count", "annot"))
filenames.list[["C3_16s"]] <- setNames(object = c("NGS5-16Schord_S1_L001_ali_assi_uniq_c10_55-775_clean_HS.txt", "NGS5-16Schord_S1_L001_ali_assi_uniq_c10_55-775_clean_HS-ex_sp.txt")
                                       , nm = c("count", "annot"))
filenames.list[["C3_COI"]] <- setNames(object = c("NGS5-C01sal_S2_L001_ali_assi_uniq_c10_55-775_clean_HS.txt", "NGS5-C01sal_S2_L001_ali_assi_uniq_c10_55-775_clean_HS_annot-ex_sp.txt")
                                       , nm = c("count", "annot"))
filenames.list[["SOG_val"]] <- setNames(object = c("all_files_ali_assi_uniq_c10_55-75_clean_HS.txt", "all_files_ali_assi_uniq_c10_55-75_clean_HS_annot-ex_sp.txt")
                                        , nm = c("count", "annot"))
filenames.list[["SOG_16s"]] <- setNames(object = c("NGS4-16Schord_S1_L001_ali_assi_uniq_c10_55-775_clean_HS.txt", "NGS4-16Schord_S1_L001_ali_assi_uniq_c10_55-775_clean_HS-1-ex_sp.txt")
                                        , nm = c("count", "annot"))
filenames.list[["SOG_COI"]] <- setNames(object = c("all_files_ali_assi_uniq_c10_55-775_clean_HS.txt", "all_files_ali_assi_uniq_c10_55-775_clean_HS_hits-ex_sp.txt")
                                        , nm = c("count", "annot"))


#### 1.0 Import input data and merge #####
paste("You are analyzing ", datatype, sep = "")
counts <- read.delim2(paste("04_samples/", filenames.list[[datatype]][1], sep = "")) 
annot <- read.delim2(paste("05_annotated/", filenames.list[[datatype]][2], sep = ""), header = F
                     , col.names = c("id","taxon"))
# head(counts)
# head(annot)
names(counts) 
names(annot)

# Sort data
counts.sorted <- counts[order(counts$id),]
annot.sorted <- annot[order(annot$id),]

# Merge
data <- merge(x = annot.sorted, y = counts.sorted, by = "id")
names(data)

# Count up the reads coming out of MEGAN
sum(data$count)

##### 1.1 Limit to amplicon size #####
# Investigate which species that you will remove with a particular amplicon size filter
losing.species <- sort(unique(data[which(data$seq_length > 250), "taxon"]))
losing.species
keep.species <- sort(unique(data[which(data$seq_length <= 250), "taxon"]))
keep.species

# what species in the losing species are also present in the keep species?
losing.species[losing.species %in% keep.species]

# Calculate how many reads are lost when filtering by size here
# data[which(data$seq_length > 250), ]
remove.sum <- sum(data[which(data$seq_length > 250), "count"]) # removing this
keep.sum <- sum(data[which(data$seq_length <= 250), "count"]) # keeping this
remove.sum / (remove.sum + keep.sum) * 100 # percent being lost by this filter

# Specific species
target.species <- "Squalus acanthias"
target.sp.removed.sum <- sum(data[data$taxon==target.species & data$seq_length > 250, "count"])
target.sp.retained.sum <- sum(data[data$taxon==target.species & data$seq_length <= 250, "count"]) 
total.sum <- sum(data[data$taxon==target.species, "count"])
target.sp.removed.sum / total.sum * 100

# Finally, limit to the rows with amplicons <= 250 bp
data2 <- data[data$seq_length <= 250, ]
sum(data2$count)
sum(data$count)
data <- data2
str(data)

###### 1.2 Reduce columns within data ####
data.df <- as.data.frame(data[, grepl( "sample\\.|taxon", names( data ))]) # keeps 'sample.' or 'taxon'
head(data.df)

# View species that are present in dataset
unique(data.df$taxon)


#### 1.3 Set location information ####
locations <- list()
locations[["locations.C3"]] <- c("IleQuarry", "Charlott", "LouisbNS", "TerraNova","RigolNL","RamahNL"
                                 , "PondInlet" , "ErebusNu", "StRochNu", "BathhurNu", "PearceNT", "NomeAK"
                                 , "HaidaGwaiiBC", "KutzeBC",  "ExtCont", "NTC")
locations[["locations.SOG"]] <- colnames(data.df)[2:length(colnames(data.df))]

location.type <- paste("locations.", sub("_\\S*", "", datatype), sep = "")

sample.locations <- locations[[location.type]]
sample.locations


##### 1.4  Explore unassigned or unannotated data #####
head(data.df)
table.filename <- paste("05_annotated/", datatype, "_unassigned_unknown_counts.csv", sep = "")
no.hits <- colSums(data.df[data.df$taxon == "No hits" , 2:length(colnames(data.df))])
not.assigned <- colSums(data.df[data.df$taxon == "Not assigned" , 2:length(colnames(data.df))])
all.hits <- colSums(data.df[, 2:length(colnames(data.df))])
percent.no.hits <- no.hits/all.hits * 100
percent.not.assigned <- not.assigned/all.hits * 100

options(scipen= 999999, digits=2)
unannot.df <- rbind(all.hits, no.hits, percent.no.hits, not.assigned, percent.not.assigned)
unannot.df <- round(x = unannot.df, digits = 2)

colnames(unannot.df) <- sample.locations

# write.csv(x = unannot.df, file = table.filename)


# Set species to remove (e.g. humans)
remove.from.all <- c("Not assigned", "No hits")
species.remove <- list()
species.remove[["C3_16s"]] <- c("Homininae", "Homo sapiens", remove.from.all)
species.remove[["C3_COI"]] <- c("NA", remove.from.all)
species.remove[["C3_val"]] <- c("Homo sapiens", "Homininae", remove.from.all)
species.remove[["SOG_val"]] <- c("Homo sapiens", remove.from.all)
species.remove[["SOG_16s"]] <- c("Homo sapiens", "Homininae", remove.from.all)
species.remove[["SOG_COI"]] <- c("Homo sapiens", "Homininae", remove.from.all)
species.remove <- species.remove[[datatype]] # Use datatype for removing species
species.remove

# Remove species from data.df
dim(data.df)
data.df <- data.df[ ! data.df$taxon %in% species.remove, ]
dim(data.df) # see how the number of taxa is reduced




#### 2. Get proportional and count data by taxon per site ####
# Set nulls
sample.tot <- NULL ; sample <- NULL; result <- NULL; result.prop <- NULL ; count.list <- NULL
prop.list <- list(); agg.counts.list <- list()

# Loop to get count by species by site (count.list) and proportion of species by site (prop.list)
for(col in 2:ncol(data.df)) {
  
  # name of the sample this iteration
  sample <- names(data.df[col]) 
  # total number reads in this sample
  sample.tot <- sum(data.df[,col])
  
  # Add 0.0001 to avoid 0 value if sample.tot is 0 (e.g. for controls)
  if(sample.tot==0){
    sample.tot <- 0.00001
  }
  
  # Per sample, aggregate counts by taxon
  result <- aggregate(x = data.df[,col], by = list(data.df$taxon), FUN = sum, na.rm = T)
  
  # Make result proportional by dividing by the amount for that species (2nd column in 'result') 
  result.prop <- (result[,2]/sample.tot)*100 
  
  # Save count values into a list
  count.list[[sample]] <- setNames(result[,2], result[,1]) # pull the value into count.list with names as the first column of result
  
  # Save proportions into a list
  prop.list[[sample]] <- setNames(result.prop, result[,1])
}
  
str(prop.list)
str(count.list)

#### 3. Pull out relevant information from proportional and count data #####
# Proportion data
prop.df <- NULL 

for(i in 1:length(prop.list)){ 
  prop.df <- cbind(prop.df, prop.list[[i]])}

head(prop.df)

# Count data
counts.df <- NULL

for(i in 1:length(count.list)){ 
  counts.df <- cbind(counts.df, count.list[[i]])}

head(counts.df)

# Incorporate location names
site.names <- sample.locations[1:length(prop.df[1,])] # Assumes is in same order for the two major types individually (SOG or C3)

colnames(prop.df) <- site.names # name prop.df
colnames(counts.df) <- site.names # name counts.df
head(prop.df)
head(counts.df)

# Set filenames for saving out
count.output.csv.filename <- paste("05_annotated/", datatype, "_count_by_taxa.csv", sep = "")
prop.output.csv.filename <- paste("05_annotated/", datatype, "_prop_by_taxa.csv", sep = "")
# write.csv(x = counts.df, file = count.output.csv.filename)
# write.csv(x = prop.df, file = prop.output.csv.filename)

# Find total numbers of reads mapping
colnames(counts.df)
sample.reads <- colSums(x = counts.df[, c(1:ncol(counts.df))])


# Filter counts table to remove very low values
min.count <- 10
counts.filtered.df <- counts.df

head(counts.filtered.df)
counts.filtered.df[which(counts.filtered.df < min.count)]

counts.filtered.df[which(counts.filtered.df < min.count)] <- 0
head(counts.filtered.df)

counts.filtered.filename <- paste("05_annotated/", datatype, "_count_by_taxa_filt_at_", min.count, ".csv", sep = "")
# write.csv(x = counts.filtered.df, file = counts.filtered.filename)

##### 4.0 Prepare plotting (colors) ####
# Prepare palette
#display.brewer.all() # see color options
cols <- brewer.pal(n = 9, name = "Set1")
cols2 <- brewer.pal(n = 8, name = "Set2")
cols3 <- brewer.pal(n = 10, name = "Set3")
cols4 <- brewer.pal(n = 8, name = "Pastel2")
cols5 <- brewer.pal(n = 9, name = "Pastel1")
cols6 <- brewer.pal(n = 11, name = "BrBG")
cols7 <- brewer.pal(n = 10, name = "Paired")
cols8 <- brewer.pal(n = 11, name = "Spectral")
cols9 <- brewer.pal(n = 9, name = "YlOrRd")
cols10 <- brewer.pal(n = 9, name = "YlGnBu")
cols11 <- brewer.pal(n = 9, name = "YlGn")
cols12 <- brewer.pal(n = 9, name = "RdPu")
cols13 <- brewer.pal(n = 9, name = "Purples")
cols14 <- brewer.pal(n = 9, name = "PuRd")
cols15 <- brewer.pal(n = 9, name = "Greys")
cols16 <- brewer.pal(n = 11, name = "RdGy")

palette <- c(cols,cols2,cols3,cols4,cols5,cols6,cols7,cols8,cols9,cols10,cols11,cols12,cols13,cols14,cols15,cols16)
length(palette)

# Randomly select from palette
set.seed(123)
index <- sample(1:nrow(counts.df))
index

# In case need numerous sets
palette.numerous<- rep(x = palette, times = 4)

# Set up in case the number is greater than a single palette
if(length(index) > length(palette)){
  this.palette <- palette.numerous[index]
} else {
  this.palette <- palette[index]
}


#### 4.1 Drop Mock Column ####
# This will remove the mock column for purposes of plotting as it overwhelms all of the data (due to too much sequencing for this)
if("sample.Mock" %in% colnames(counts.df) == T){
  counts.df <- counts.df[,-(which(colnames(counts.df)=="sample.Mock"))]
  prop.df <- prop.df[,-(which(colnames(prop.df)=="sample.Mock"))]
  site.names <- site.names[-(which(site.names == "sample.Mock"))]
  sample.reads <- sample.reads[-(which(names(sample.reads)=="sample.Mock"))]
}


##### 4.2 Create Legend ####
# Prepare legend size 
legend.cex <- c(1, 1, 1, 1, 1, 1) ; names(legend.cex) <- c("C3_16s","C3_COI", "SOG_16s", "C3_val", "SOG_val", "SOG_COI")

# Create dataframe with the taxon and the color
color.index <- cbind(rownames(prop.df), this.palette)
colnames(color.index) <- c("taxon","color")
head(color.index)
color.index.df <- as.data.frame(color.index)
# note this will not work until you move the color codes to character

# Identify which taxa are high proportion in any one sample
min.proport <- 5

# Set null
high.presence.taxa <- NULL

for(i in 1:nrow(prop.df)){
  high.presence.taxa.add <- if(max(prop.df[i,]) > min.proport) { print(rownames(prop.df)[i])}
  high.presence.taxa <- c(high.presence.taxa, high.presence.taxa.add)
}

high.presence.taxa

# Select the rows of the color index for only those high.presence taxa
legend.info <- color.index.df[color.index.df$taxon %in% high.presence.taxa, ]


#### 5. Plot  ####
filename <- paste("06_output_figures/", datatype, "_read_count_and_prop_by_loc.pdf", sep = "")

# if want to work interactively, comment out the following line
pdf(file = filename, width = 10, height = 8)
par(mfrow=c(3,1), mar= c(2.5,4.5,2,1) + 0.2, mgp = c(3.75, 0.75, 0))

# Plot count data
position.info <- barplot(as.matrix(counts.df)
                         , col = this.palette, las = 2, xaxt = "n"
                         , xlim = c(0, ncol(prop.df)+4)
                         , ylab = "Reads")
# axis(side = 1, at = position.info, labels = sample.locations, las = 3, cex.axis = 0.9)

# Plot proportion data
position.info <- barplot(as.matrix(prop.df), col = this.palette
        , xlim = c(0, ncol(prop.df)+4)
        , las = 1
        , ylab = "Proportion (%)"
        , xaxt = "n")

axis(side = 1, at = position.info, 
     labels = site.names, las = 3
     )

# Add information about read counts per sample
mtext(x = position.info, text = sample.reads
      , side=3, at = position.info, cex = 0.7)


# Plot Legend, first make blank second plot
plot(1, type = "n", axes = F, xlab = "", ylab = "")

# fix legend info to character text
legend(x = "center", y = "center", legend = legend.info$taxon
        , fill = as.character(legend.info$color), cex = legend.cex[datatype]
        , ncol = 4)

dev.off()
#
# Save out as 10 x 8 in portrait

