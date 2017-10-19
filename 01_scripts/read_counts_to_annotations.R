# Use this script to connect the read counts to the annotations
setwd("~/Documents/03_eDNA/eDNA_metabarcoding_SOG_val")

# Install Packages
#install.packages("ggplot2")
library("ggplot2")




# Import data
counts <- read.delim2("04_samples/NGSLib1_cleanHS.txt") 
annot <- read.delim2("05_annotated/NGSLib1_cleanHS_hits-ex_sp.txt", header = F
                     , col.names = c("id","species"))

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
data.df <- as.data.frame(data[,c(2,4,5:10)])
head(data.df)


# Make data into proportional data
sample.tot <- NULL ; sample <- NULL
result <- NULL
result.prop <- NULL
result.list <- list()

for(col in 3:8) {
  sample <- names(data.df[col])
  sample.tot <- sum(data.df[,col])
  result <- aggregate(x = data.df[,col], by = list(data.df$species), FUN = sum, na.rm = T)
  
  # divide the sum for that species (2nd column) by the total to make proportional
  result.prop <- (result[,2]/sample.tot)*100  
  result.list[[sample]] <- result.prop
}
  
str(result.list)

# Obtain data from list
test <- NULL 
for(i in 1:6){ 
  test <- cbind(test, result.list[[i]])}

colnames(test) <- gsub(pattern = "sample.", replacement = "S_", x = names(result.list))
rownames(test) <- result[,1] # this is too hacky
test



# Find total numbers of reads mapping
colnames(data.df)
sample.reads <- colSums(x = data.df[,c(3:ncol(data.df))])




barplot(as.matrix(test))

#install.packages("RColorBrewer")
library("RColorBrewer")
# see color options
#display.brewer.all()
# HACKY!
cols <- brewer.pal(n = 10, name = "Set1")
cols2 <- brewer.pal(n = 10, name = "Set2")
cols3 <- brewer.pal(n = 10, name = "Set3")
length(c(cols, cols2, cols3))
pallette <- c(cols,cols2,cols3)


# PLOT
# par(mfrow=c(1,1), mar= c(4,3,0.5,1) + 0.2, mgp = c(2,0.75,0))
position.info <- barplot(as.matrix(test), col = pallette
        , xlim = c(0, ncol(test)+4)
        #, legend.text = rownames(test)
        , las = 1
        , cex.names = 0.9
        , cex.axis = 0.9
        , ylab = "Proportion (%)")

legend("bottomright", legend = rownames(test)
       , fill = pallette, cex = 0.7)

test

text(x = position.info
       #1:ncol(test)
       , y = 104, labels = sample.reads, cex = 0.7)

#











## GGPLOT


library(reshape)
test2 <- melt(cbind(test, ind = rownames(test)), id.vars = c('ind'))
#




dat <- read.table(text = "    ONE TWO THREE
1   23  234 324
2   34  534 12
3   56  324 124
4   34  234 124
5   123 534 654",sep = "",header = TRUE)

#Add an id variable for the filled regions
library(reshape)
datm <- melt(cbind(dat, ind = rownames(dat)), id.vars = c('ind'))

library(scales)
ggplot(datm,aes(x = variable, y = value,fill = ind)) + 
  geom_bar(position = "fill",stat = "identity") +
  # or:
  # geom_bar(position = position_fill(), stat = "identity") 
  scale_y_continuous(labels = percent_format())
#


ggplot(data = data, aes(x = data$sample.2 ))


