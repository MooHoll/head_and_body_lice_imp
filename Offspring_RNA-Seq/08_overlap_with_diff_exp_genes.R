###-------------------------------------------------------------------
# Overlap of imprinted/ecotype genes with diff exp genes between pure lines
###-------------------------------------------------------------------

setwd("~/Dropbox/Leicester_postdoc/Projects/lice/RNA_seq_pure")

library(readr)
library(ggplot2)

imp_genes_output <- read_delim("~/Dropbox/Leicester_postdoc/Projects/lice/imp_genes_output.txt", 
                               delim = "\t", escape_double = FALSE, 
                               trim_ws = TRUE)
all_gene_expression_data <- read_delim("counts/Pure_data/all_gene_expression_data.txt", 
                                       delim = "\t", escape_double = FALSE, 
                                       trim_ws = TRUE)

###-------------------------------------------------------------------

# Take only the interesting maternal/ecotype specific genes
sig_imp_genes <- imp_genes_output[imp_genes_output$actually_sig=="yes",]
sig_imp_genes <- sig_imp_genes[!is.na(sig_imp_genes$hb),] # remove that one weird gene again
table(sig_imp_genes$lineage) # This is right

with_diff_exp <- merge(sig_imp_genes, all_gene_expression_data, by="gene_id") # all 73 are there

# Ok, almost none of them are differentially expressed between adult males of the pure lines
# One: gene-Phum_PHUM549420 which shows head lice allele expression bias in the hybrids of around 73-78%
# but shows higher expression in pure body lice compared to head lice with an SPM of 0.92
# gene name: "Programmed cell death protein, putative".

###-------------------------------------------------------------------
# Plot the difference in expression for the imp genes anyway



