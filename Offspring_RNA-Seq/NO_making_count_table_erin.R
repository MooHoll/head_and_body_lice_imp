###-------------------------------------------------------------------
# Making input file for the imprinting model
###-------------------------------------------------------------------

setwd("~/Downloads")

library(readr)
library(data.table)
library(dplyr)
library(doBy)

###-------------------------------------------------------------------
# Read in sample count files for maternal reads

# Make a list of the files you want
file.list = list.files(("./"),pattern="*f_counts.txt")

# Create a function to read them in so you can put all files in a list
read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = F, trim_ws = T)
}

# Actually read in the files
samples <- lapply(file.list, read_file1)

# Create a vector with the file names 
sample_names <- tools::file_path_sans_ext(file.list)

# Take only the first part of the file name for the sample name
sample_names <- lapply(X = sample_names, FUN = function(t) gsub(pattern = "_f_counts", 
                                               replacement = "", x = t, fixed = TRUE))

# Name the files in the list with their sample name
names(samples) <- sample_names

# Add a column for the sample name and name the existing columns
# Also get rid of the extra annotation information that isn't a gene name
# Also for maternal reads if they're zero we need to remove them completely (see methods in Galbraith et al. 2015)
for(i in seq_along(samples)){
  colnames(samples[[i]]) <- c("gene_id", "f_count")
  samples[[i]]$sample_name <- names(samples[i])
  samples[[i]] <- samples[[i]][!samples[[i]]$gene_id %like% "__",]
  samples[[i]] <- samples[[i]][!samples[[i]]$f_count == 0,]
}

# Put all the maternal counts into one table
maternal_all <- as.data.frame(bind_rows(samples))

###-------------------------------------------------------------------
# Do the same for the paternal reads

file.list = list.files(("./"),pattern="*m_counts.txt")
read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = F, trim_ws = T)
}
samples <- lapply(file.list, read_file1)
sample_names <- tools::file_path_sans_ext(file.list)
sample_names <- lapply(X = sample_names, FUN = function(t) gsub(pattern = "_m_counts", 
                                                               replacement = "", x = t, fixed = TRUE))
names(samples) <- sample_names

# The only difference here is if males have 0 we convert to 1 to avoid complete separation in the model
for(i in seq_along(samples)){
  colnames(samples[[i]]) <- c("gene_id", "m_count")
  samples[[i]]$sample_name <- names(samples[i])
  samples[[i]] <- samples[[i]][!samples[[i]]$gene_id %like% "__",]
  samples[[i]]$m_count[samples[[i]]$m_count==0] <- 1
}
paternal_all <- as.data.frame(bind_rows(samples))

###-------------------------------------------------------------------
# Put it all together

# Merge both dataframes
all_data <- merge(maternal_all, paternal_all, by=c("gene_id","sample_name"))

# Make sure we only keep genes which have data for all samples
all_data <- all_data[complete.cases(all_data),]

# Add a total count row
all_data$total <- all_data$m_count+all_data$f_count

# Remove genes with less total coverage than 5 for sample
all_data <- all_data[all_data$total>5,] 

# See how many genes we're actually testing in the model
length(unique(all_data$gene_id))


# TO DO!!!! Add in a column with the direction of the cross. If you have a file
# which has the sample name and direction of cross as two columns
# you can read it in and do a merge command to add the cross direction
# to the all_data file

# As I don't have this file I'm going to make it up for the below
sample_name <- c("SRR2773794","SRR2773795","SRR2773796","SRR2773797")
cross_direction <- c("VG", "GV","VG","GV")
cross_info <- data.frame(sample_name, cross_direction)

all_data_with_cross <- merge(all_data, cross_info)

###-------------------------------------------------------------------
# Only keep genes present in both cross directions

# Make a file with genes per cross, irrespective of sample number (think about this later)
genes_with_cross <- all_data_with_cross[,c(2,6)]
genes_with_cross <- genes_with_cross[!duplicated(genes_with_cross),]

# See now many crosses each gene is in (either one or two)
genes_in_both_crosses <- summaryBy(cross_direction ~ gene_id, data=genes_with_cross, FUN=length)

# Make a list of genes which are in both crosses
genes_in_both_crosses <- genes_in_both_crosses$gene_id[genes_in_both_crosses$cross_direction==2]
length(genes_in_both_crosses)

# Only keep genes present in both crosses
all_data_with_cross <- all_data_with_cross[all_data_with_cross$gene_id %in% genes_in_both_crosses,]

###-------------------------------------------------------------------
# Write out a dataframe for the model
write.table(all_data_with_cross, file ="all_data_for_imp_model.txt", sep = "\t",
            quote = F, col.names = T, row.names = F)




