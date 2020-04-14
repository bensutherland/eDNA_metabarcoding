# Rarefy read counts to plot taxa counts based on sample with highest taxa number
# rm(list=ls())
setwd("~/03_eDNA/01_manuscripts_and_projects/2020_C3_pilot_and_SOG/04_results_and_analyses/read_rarefaction")

# Select your marker (i.e., 16s, COI, val)
# marker <- "16s"
# marker <- "COI"
# marker <- "val"


#### 00. Set variables ####
plot_type <- "all" # indiv or all
input.PATH <- "C:/00_GitHub/eDNA_metabarcoding/07_results/"
sample_type <- "random" # "random" or "abund"
drop_number <- 200 # 50 reads will be dropped from the total


# # Plot a single pane or three pane
# if(plot_type=="indiv"){
#   
#   par(mfrow = c(1,1))
# }else if(plot_type=="all"){
#   
#   par(mfrow = c(3,1))
#   
# }

pdf(file = "rarefaction.pdf", width = 5, height = 8)
par(mfrow = c(3,1))

# Loop across the three markers
markers <- c("16s", "COI", "val")
moi <- NULL

for(m in 1:length(markers)){
  
  # Select the marker of interest
  moi <- markers[m]


    # Read in data
    input.FN <- list.files(path = input.PATH, pattern = moi)
    print(paste0("Reading in ", input.FN))
    
    data <- read.csv(file = paste0(input.PATH, input.FN), row.names = 1)
    head(data)
    
    
    # Identify sample with the most abundant taxa
    which(colSums(x = data, na.rm = T)==max(colSums(x = data, na.rm = T))) # most reads
    most_abund_sample <- names(sort(x = colSums(data != 0), decreasing = T))[1] # most taxa
    
    
    # Either select the most abundant sample or a random
    if(sample_type=="random"){
      
      most_abund_sample <- names(data)[sample(x = seq(1:length(names(data))), size = 1)]
      print(paste0("using a random sample: ", most_abund_sample))
      
    }else if(sample_type=="abund"){
      
      print(paste0("using the most abundant sample: ", most_abund_sample))
      
    }
    
    
    # Remove specific unwanted species
    sp_to_rem <- c("cellular organisms", "No hits", "Not assigned") # specific to 16S currently
    `%notin%` <- Negate(`%in%`)
    data <- data[rownames(data) %notin% sp_to_rem, ]
    
    
    # Reduce the dataset to only the chosen site
    data_to_rarefy.df <- as.data.frame(data[,most_abund_sample])
    colnames(data_to_rarefy.df) <- most_abund_sample
    head(data_to_rarefy.df)
    str(data_to_rarefy.df)
    sum(data_to_rarefy.df)
    
    
    #### Rarefy reads #####
    # Set up list to capture results
    taxa.list <- list(); total_reads.list <- list()
    iterations <- NULL
    
    # reduce to non-zero
    sample_taxa.dat <- data_to_rarefy.df[data_to_rarefy.df!=0, most_abund_sample]
    head(sample_taxa.dat)
    sum(sample_taxa.dat)
    str(sample_taxa.dat)
    # sample_taxa.dat is now an integer vector
    
    # determine number of iterations needed
    iterations <- sum(sample_taxa.dat) / drop_number
    iterations <- round(x = iterations, digits = 0)
    
    
    # Loop to remove values iteratively
    for(j in 1:iterations){
      
      # Subtract n reads from random taxa counts
      for(i in 1:drop_number){
        
        reduce_record <-  sample(x = seq(1:length(sample_taxa.dat)), size = 1)
        
        # Subtract 1 from the selected record
        sample_taxa.dat[reduce_record] <- sample_taxa.dat[reduce_record] - 1
        
        # Remove zero values from the vector
        sample_taxa.dat <- sample_taxa.dat[sample_taxa.dat!=0]
        
      }
      
      print(paste0("total reads left: ", sum(sample_taxa.dat)))
      
      # Collect data for this iteration
      taxa.list[[j]] <- length(sample_taxa.dat)
      total_reads.list[[j]] <- sum(sample_taxa.dat)
      iterations_count <- seq(from = 1, to = iterations, by = 1)
      
    }
    
    length(taxa.list)
    length(total_reads.list)
    length(iterations_count)
    
    #### Plot Results ####
    
    # Ascending order:
    taxa_count <- unlist(x = taxa.list)
    read_count <- unlist(x = total_reads.list)
    # and use iterations_count
    
    output.df <- as.data.frame(
                        cbind(taxa_count, read_count, iterations_count)
                        , stringsAsFactors = F) 
    head(output.df)
    
    
    plot(x = output.df$iterations_count, y = output.df$taxa_count, main = paste0("rarefy ", moi, "- ", most_abund_sample)
         , xlab = "Number Reads (x1000)", ylab = "Number Taxa"
         , type = "l"
         , las = 1
         , xaxt = "n"
         )
    total_tenth <- round(length(iterations_count) / 10, digits = 0)
    axis(side = 1, at = output.df$iterations_count[seq(1, nrow(output.df), total_tenth)]
         , labels = round(output.df$read_count[seq(1, nrow(output.df), total_tenth)] / 1000, digits = 0), las = 2)
    
    }

dev.off()
