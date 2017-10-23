# Use this script to connect the read counts to the annotations
#rm(list=ls())


#setwd("~/Documents/03_eDNA/eDNA_metabarcoding_SOG_val")
#setwd("~/Documents/03_eDNA/eDNA_metabarcoding_C3_16s")
setwd("~/Documents/03_eDNA/eDNA_metabarcoding_C3_COI")

# Install Packages
#install.packages("RColorBrewer")
library("RColorBrewer")

# Choose datatype
datatype <- "C3_16s"
#datatype <- "C3_COI"

## Set filenames for loading different datasets
filenames.list <- list()
filenames.list[["C3_16s"]] <- c("NGS5_C3_16S_cleanHS.txt", "NGS5_C3_16S_cleanHS_hits-1-ex_species.txt")
names(filenames.list[["C3_16s"]]) <- c("count", "annot")
filenames.list[["C3_COI"]] <- c("NGS_C3_cleanHS.txt", "NGS_C3_COI_cleanHS_hits-ex_species.txt")
names(filenames.list[["C3_COI"]]) <- c("count", "annot")

filenames.list


# Import data
paste("You are analyzing ", datatype, sep = "")
counts <- read.delim2(paste("04_samples/", filenames.list[[datatype]][1], sep = "")) 
annot <- read.delim2(paste("05_annotated/", filenames.list[[datatype]][2], sep = ""), header = F
                     , col.names = c("id","taxon"))

head(counts)
head(annot)

names(counts)
names(annot)

# sort
counts.sorted <- counts[order(counts$id),]
annot.sorted <- annot[order(annot$id),]

# merge
data <- merge(x = annot.sorted, y = counts.sorted, by = "id")
names(data)
str(data)


###### Quantification #####
# subset only the required columns
data.df <- as.data.frame(data[, grepl( "sample\\.|taxon", names( data ))]) # more adaptive method
head(data.df)

# Species remove (#TODO)

# Make data proportional data
sample.tot <- NULL ; sample <- NULL; result <- NULL; result.prop <- NULL ; result.count.list <- NULL
result.list <- list(); agg.counts.list <- list()

for(col in 2:ncol(data.df)) {
  sample <- names(data.df[col]) # name of the sample this iteration
  sample.tot <- sum(data.df[,col]) # total number reads in this sample
  
  # Add up all of the counts for this sample by taxon
  result <- aggregate(x = data.df[,col], by = list(data.df$taxon), FUN = sum, na.rm = T)
  
  # divide the sum for that species (2nd column in 'result') 
  # by the total reads for that sample to make proportional
  result.prop <- (result[,2]/sample.tot)*100 
  
  # save out counts into a list
  result.count.list[[sample]] <- result[,2]
  # save out proportions into a list
  result.list[[sample]] <- result.prop
}
  
str(result.list)
str(result.count.list)

#### Obtain relevant information into a dataframe from the list ####
# for proportions
prop.df <- NULL 
for(i in 1:length(result.list)){ 
  prop.df <- cbind(prop.df, result.list[[i]])}

colnames(prop.df) <- gsub(pattern = "sample.", replacement = "S_", x = names(result.list))
rownames(prop.df) <- result[,1] # a bit too hacky
head(prop.df)
# write.csv(x = prop.df, file = "05_annotated/proportions_by_taxa.csv")

# for counts
counts.df <- NULL
for(i in 1:length(result.count.list)){ 
  counts.df <- cbind(counts.df, result.count.list[[i]])}

colnames(counts.df) <- gsub(pattern = "sample.", replacement = "S_", x = names(result.list))
rownames(counts.df) <- result[,1] # a bit hacky
head(counts.df)
# write.csv(x = counts.df, file = "05_annotated/counts_by_taxa.csv")
##TODO## MAKE THAT MORE ADAPTABLE BY USING PASTE

# Find total numbers of reads mapping
colnames(data.df)
sample.reads <- colSums(x = data.df[, c(2:ncol(data.df))])


##### Plot proportion data ####

# see color options
#display.brewer.all()
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

# Here we repeat colors, because of so many species (repeats will probably not be viewed as only a subset will be shown in the legend in the end)
palette.numerous<- rep(x = palette, times = 4)


### Connect locations to plot ####
locations.list <- list()
locations.list[["C3_16s"]] <- c("IleQuarry", "Charlott", "LouisbNS", "TerraNova","RigolNL","RamahNL"
                                             , "PondInlet" , "ErebusNu", "StRochNu", "BathhurNu", "PearceNT", "NomeAK"
                                             , "HaidaGwaiiBC", "KutzeBC",  "ExtCont", "NTC")
locations.list[["C3_COI"]] <- c("IleQuarry", "Charlott", "LouisbNS", "TerraNova","RigolNL","RamahNL"
                                , "PondInlet" , "ErebusNu", "StRochNu", "BathhurNu", "PearceNT", "NomeAK"
                                , "HaidaGwaiiBC", "KutzeBC",  "ExtCont")

# Select locations based on datatype
locations <- locations.list[[datatype]]

# Make index for sample names and locations
sample.locations <- as.data.frame(cbind(colnames(prop.df), locations))


# Minor adjust needed for cex in legend for the 16s or COI
legend.cex <- c(1.2, 0.9) ; names(legend.cex) <- c("C3_16s","C3_COI")


# PLOT
par(mfrow=c(2,1), mar= c(4,3,3,1) + 0.2, mgp = c(2,0.75,0))
position.info <- barplot(as.matrix(prop.df), col = palette.numerous[1:nrow(prop.df)]
        , xlim = c(0, ncol(prop.df)+4)
        , las = 1
        , cex.names = 0.9
        , cex.axis = 0.9
        , ylab = "Proportion (%)"
        , xaxt = "n")

axis(side = 1, at = position.info, 
     labels = sample.locations$locations, las = 3
     , cex.axis = 0.9)

# Add information about read counts per sample
# text(x = position.info
#      , y = 140, labels = sample.reads, cex = 0.7)
mtext(x = position.info, text = sample.reads
      , side=3, at = position.info, cex = 0.7)


# Create a legend index to only create legend for top presence
color.index <- cbind(rownames(prop.df), palette.numerous[1:nrow(prop.df)])
colnames(color.index) <- c("taxon","color")
head(color.index)
color.index.df <- as.data.frame(color.index)
# note this will not work until you move the color codes to character

# make a legend with only those with greater than 5% contained
head(prop.df)

# Identify which taxa are high proportion in any one sample
min.proport <- 5

high.presence.taxa <- NULL
for(i in 1:nrow(prop.df)){
  high.presence.taxa.add <- if(max(prop.df[i,]) > min.proport) { print(rownames(prop.df)[i])}
  high.presence.taxa <- c(high.presence.taxa, high.presence.taxa.add)
}

high.presence.taxa

legend.info <- color.index.df[color.index.df$taxon %in% high.presence.taxa, ]


# blank second plot
plot(1, type = "n", axes = F, xlab = "", ylab = "")

# fix legend info to character text
legend(x = "center", y = "center", legend = legend.info$taxon
        , fill = as.character(legend.info$color), cex = legend.cex[datatype]
        , ncol = 3)

#
# Save out as 10 x 8 in portrait
