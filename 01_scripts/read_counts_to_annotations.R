# Use this script to connect the read counts to the annotations
#setwd("~/Documents/03_eDNA/eDNA_metabarcoding_SOG_val")
setwd("~/Documents/03_eDNA/eDNA_metabarcoding_C3_16s/")

# Install Packages
#install.packages("RColorBrewer")
library("RColorBrewer")


# Set files
count.file.name <- "NGS5_C3_16S_cleanHS.txt"
annot.file.name.class <- "NGS5_C3_16S_cleanHS_hits-1-ex_order.txt"
annot.file.name.order <- "NGS5_C3_16S_cleanHS_hits-1-ex_order.txt"
annot.file.name.species <- "NGS5_C3_16S_cleanHS_hits-1-ex_species.txt" 

annot.file.name <- annot.file.name.class

# Import data
counts <- read.delim2(paste("04_samples/", count.file.name, sep = "")) 
annot <- read.delim2(paste("05_annotated/", annot.file.name, sep = ""), header = F
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
#data.df <- as.data.frame(data[,c(2,4,5:10)])
head(data.df)

# Make data proportional data
sample.tot <- NULL ; sample <- NULL
result <- NULL
result.prop <- NULL
result.list <- list()

for(col in 2:ncol(data.df)) {
  sample <- names(data.df[col])
  sample.tot <- sum(data.df[,col])
  result <- aggregate(x = data.df[,col], by = list(data.df$taxon), FUN = sum, na.rm = T)
  
  # divide the sum for that species (2nd column) by the total to make proportional
  result.prop <- (result[,2]/sample.tot)*100  
  result.list[[sample]] <- result.prop
}
  
str(result.list)


# Bring data out of list into a dataframe
prop.df <- NULL 
for(i in 1:length(result.list)){ 
  prop.df <- cbind(prop.df, result.list[[i]])}

colnames(prop.df) <- gsub(pattern = "sample.", replacement = "S_", x = names(result.list))
rownames(prop.df) <- result[,1] # a bit too hacky
head(prop.df)


# Find total numbers of reads mapping
colnames(data.df)
sample.reads <- colSums(x = data.df[, c(2:ncol(data.df))])


##### PLOT ####
# see color options
#display.brewer.all()

# option 1
#palette <- colorRampPalette(c("blue","red"))(nrow(prop.df))

# option 2
cols <- brewer.pal(n = 10, name = "Set1")
cols2 <- brewer.pal(n = 10, name = "Set2")
cols3 <- brewer.pal(n = 10, name = "Set3")
cols4 <- brewer.pal(n = 8, name = "Pastel2")
cols5 <- brewer.pal(n = 9, name = "Pastel1")
cols6 <- brewer.pal(n = 10, name = "Paired")
cols7 <- brewer.pal(n = 11, name = "BrBG")
cols8 <- brewer.pal(n = 11, name = "Spectral")
palette <- c(cols,cols2,cols3,cols4,cols5,cols6,cols7,cols8)
length(palette)


### Connect locations to plot ####
locations <- c("IleQuarry", "Charlott", "LouisbNS", "TerraNova","RamahNL","RigolNL"
               , "PondInlet" , "ErebusNu", "StRochNu", "BathhurNu", "PearceNT", "NomeAK"
               , "KutzeBC", "HaidaGwaiiBC", "ExtCont", "NTC")
sample.locations <- as.data.frame(cbind(colnames(prop.df), locations))



# PLOT
par(mfrow=c(2,1), mar= c(4,3,3,1) + 0.2, mgp = c(2,0.75,0))
position.info <- barplot(as.matrix(prop.df), col = palette[1:nrow(prop.df)]
        , xlim = c(0, ncol(prop.df)+4)
        , las = 1
        , cex.names = 0.7
        , cex.axis = 0.7
        , ylab = "Proportion (%)")

# Add information about read counts per sample
# text(x = position.info
#      , y = 140, labels = sample.reads, cex = 0.7)
mtext(x = position.info, text = sample.reads
      , side=3, at = position.info, cex = 0.7)


# Create a legend index to only create legend for top presence
color.index <- cbind(rownames(prop.df), palette[1:nrow(prop.df)])
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
        , fill = as.character(legend.info$color), cex = 0.8
        , ncol = 3)

#
