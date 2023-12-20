# ------------------------------------------------------------
# Making model input file
# ------------------------------------------------------------

setwd("~/Dropbox/Research/Leicester_postdoc/Projects/lice/snp_counts")

library(readr)
library(data.table)
library(dplyr)
library(reshape2)
library(ggplot2)
library(doBy)

# Make a list of the files you want
file.list = list.files(("./"),pattern="*genes.txt")

# Create a function to read them in so you can put all files in a list
read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = T, trim_ws = T)
}

# Actually read in the files
samples <- lapply(file.list, read_file1)

# Put all together
all_data <- as.data.frame(bind_rows(samples))

length(unique(all_data$gene_id)) #654 genes

# Keep genes in both crosses and all samples
#gene_number <- as.data.frame(table(all_data$gene_id))
#gene_number <- gene_number[gene_number$Freq==18,] #21 genes only, bloody hell

# Keep genes which have at least two replicates per condition (probably a better way than the code below)
# filter by 4 as there are two counts per replicate (hh and bb)
gene_number <- summarise(group_by(all_data,gene_id,cross),count =n())
gene_number <- spread(gene_number, key = cross, value = count)
gene_number <- gene_number[(gene_number$bh >=4 & gene_number$hb >=4),] # 509 genes

all_data <- all_data[all_data$gene_id %in% gene_number$gene_id,]

# Rename counts for maternal/paternal depending on cross
all_data$line_origin[all_data$cross =="bh" & all_data$line_origin == "bb_count"] <- "maternal" 
all_data$line_origin[all_data$cross =="bh" & all_data$line_origin == "hh_count"] <- "paternal" 

all_data$line_origin[all_data$cross =="hb" & all_data$line_origin == "hh_count"] <- "maternal" 
all_data$line_origin[all_data$cross =="hb" & all_data$line_origin == "bb_count"] <- "paternal" 

# Make back into wide format
all_data <- spread(all_data, key = line_origin, value = count.sum)

write.table(all_data, file="all_data_for_imp_model.txt", sep="\t", quote=F,
            row.names = F, col.names = T)
