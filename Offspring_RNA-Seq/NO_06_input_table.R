# ------------------------------------------------------------
# Making model input file
# ------------------------------------------------------------

setwd("~/Dropbox/Leicester_postdoc/Projects/lice/snp_counts")

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

# Add in maternal or paternal information
# Checking Andres thesis, bh = BB female with HH male and hb = HH female with BB male
for(i in seq_along(samples)){
  samples[[i]]$parent <- "maternal"
  samples[[i]]$parent[samples[[i]]$cross=="bh" & samples[[i]]$SNP_origin=="HH"] <- "paternal"
  samples[[i]]$parent[samples[[i]]$cross=="hb" & samples[[i]]$SNP_origin=="BB"] <- "paternal"
}

# Put all together
all_data <- as.data.frame(bind_rows(samples))

length(unique(all_data$gene_id)) #3980 genes

# Edit table to make it useable for parent of origin
long_all_data <- dcast(all_data, sample_name+cross+SNP_origin+chr+gene_id ~ parent, value.var = "count.sum")

# Remove lineage column as testing parent of origin
long_all_data <- long_all_data[,-3]

# Keep only genes present in both crosses
long_all_data$sample_name <- gsub('.{3}$', '', long_all_data$sample_name)

# Sort out the format
maternal <- long_all_data[,-6]
paternal <- long_all_data[,-5]

both <- merge(maternal,paternal, by=c("sample_name","cross","chr","gene_id"))
both <- both[complete.cases(both),]
length(unique(both$gene_id)) #623 genes

# Actually keep genes in both crosses and all samples
gene_number <- as.data.frame(table(both$gene_id))
gene_number <- gene_number[gene_number$Freq==9,] #384 genes only, bloody hell

both <- both[both$gene_id %in% gene_number$Var1,]

write.table(both, file="all_data_for_imp_model.txt", sep="\t", quote=F,
            row.names = F, col.names = T)

# Look
both$total <- both$maternal + both$paternal
both$proportion_mat <- both$maternal / both$total

hist(both$proportion_mat)
head(both)
for_plot <- both[,c(2,4,8)]
head(for_plot)
for_plot2 <- summaryBy(proportion_mat ~ cross + gene_id, data = for_plot, FUN = mean)

for_plot2 <- dcast(for_plot2, gene_id ~ cross, value.var = "proportion_mat.mean")

ggplot(for_plot2, aes(x = bh, y = hb))+
  geom_point()

