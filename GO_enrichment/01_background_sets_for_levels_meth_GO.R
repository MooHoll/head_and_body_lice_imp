#---------------------------------------------------
# Background gene sets for levels of meth GO analysis
#---------------------------------------------------

setwd("~/Dropbox/Leicester_postdoc/Projects/lice/GO_analysis")
library(readr)
library(tidyr)

GO_terms <- read_delim("louse_GO_annotations.txt", 
                                           delim = "\t", escape_double = FALSE, 
                                           col_names = TRUE, trim_ws = TRUE)
colnames(GO_terms) <- c("transcript","GO_id")

# Need to change gene id to actual gene id not the RNA pfft!
gtf_file <- read_delim("~/Dropbox/Leicester_postdoc/Projects/lice/genome/JCVI_LOUSE_1.0_genomic.gtf", 
                                     delim = "\t", escape_double = FALSE, 
                                     col_names = FALSE, trim_ws = TRUE)
gtf_file <- gtf_file[,9]
gtf_file <- separate(gtf_file, X9, into = c("transcript","gene_id","other"), sep = ";")
gtf_file <- gtf_file[,-3]
gtf_file$transcript <- sub(".*rna", "rna",gtf_file$transcript)
gtf_file$transcript <- sub(".1.*", ".1",gtf_file$transcript)
gtf_file$gene_id <- gsub(".{1}$", "",gtf_file$gene_id)
gtf_file$gene_id <- sub(".*gene-", "gene-",gtf_file$gene_id)

GO_for_genes <- merge(GO_terms, gtf_file, by = "transcript")
GO_for_genes <- GO_for_genes[complete.cases(GO_for_genes),]
GO_for_genes <- GO_for_genes[,-1]
length(unique(GO_for_genes$gene_id)) # 4761 only have annotations

#---------------------------------------------------
all_gene_expression_data <- read_delim("~/Dropbox/Leicester_postdoc/Projects/lice/RNA_seq_pure/all_gene_expression_data.txt", 
                                       delim = "\t", escape_double = FALSE, 
                                       trim_ws = TRUE)

all_genes <- all_genes <- as.data.frame(unique(all_gene_expression_data$gene_id))
colnames(all_genes)<-"gene_id"

all_genes_GO <- GO_for_genes[GO_for_genes$gene_id %in% all_genes$gene_id,]
length(unique(all_genes_GO$gene_id)) #4606/4761

write.table(all_genes_GO, file="rna_all_genes_pure.txt", sep="\t",quote = F, col.names = T, row.names = F)

#---------------------------------------------------
diff_exp_genes <- all_genes <- as.data.frame(unique(all_gene_expression_data$gene_id[all_gene_expression_data$diff_exp=="yes"]))
colnames(diff_exp_genes)<-"gene_id"
write.table(diff_exp_genes, file="diff_exp_genes_pure.txt", sep="\t",quote = F, col.names = T, row.names = F)

diff_exp_genes_GO <- GO_for_genes[GO_for_genes$gene_id %in% diff_exp_genes$gene_id,]
length(unique(diff_exp_genes_GO$gene_id)) #17/47

write.table(diff_exp_genes_GO, file="rna_all_diff_genes_pure.txt", sep="\t",quote = F, col.names = T, row.names = F)

#---------------------------------------------------
# Pull out other lists whilst I'm here
head(all_gene_expression_data)

BB_biased <- all_genes <- as.data.frame(unique(all_gene_expression_data$gene_id[all_gene_expression_data$category=="BB_biased"]))
colnames(BB_biased)<-"gene_id"
write.table(BB_biased, file="BB_biased_genes_pure.txt", sep="\t",quote = F, col.names = T, row.names = F)

HH_biased <- all_genes <- as.data.frame(unique(all_gene_expression_data$gene_id[all_gene_expression_data$category=="HH_biased"]))
colnames(HH_biased)<-"gene_id"
write.table(HH_biased, file="HH_biased_genes_pure.txt", sep="\t",quote = F, col.names = T, row.names = F)

BB_limited <- all_genes <- as.data.frame(unique(all_gene_expression_data$gene_id[all_gene_expression_data$category=="BB_limited"]))
colnames(BB_limited)<-"gene_id"
write.table(BB_limited, file="BB_limited_genes_pure.txt", sep="\t",quote = F, col.names = T, row.names = F)


