# Merging pvalues and frequency information 

setwd("~/Dropbox/Leicester_postdoc/Projects/lice/GO_analysis/imprinted_analysis/outputs")
library(readr)

pvals <- read_delim("maternal_expressed_limited_enriched_GOs.txt", 
                delim = "\t", escape_double = FALSE, 
                trim_ws = TRUE)

freqs <- read_delim("maternal_limited.tsv", 
            delim = "\t", escape_double = FALSE, 
            trim_ws = TRUE)

freqs <- freqs[,c(1,2,5)]
colnames(freqs)[1] <- "GOBPID"

both <- merge(pvals, freqs, by = "GOBPID")

write.table(both, file="info_for_supp.txt", sep="\t", quote = F, col.names = T, row.names = F)
