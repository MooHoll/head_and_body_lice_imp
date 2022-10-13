#------------------------------------------------------------------
# Differential Alternative Splicing
#------------------------------------------------------------------

setwd("~/Dropbox/Leicester_postdoc/Projects/lice/RNA_seq_pure/counts/Pure_data")

library(IsoformSwitchAnalyzeR)
library(readr)
library(doBy)
library(dplyr)
library(tidyr)

#------------------------------------------------------------------
# Read in files
file.list <- list.files("./", pattern = "*isoforms.results")

quant_data <- importIsoformExpression(sampleVector = file.list)

myDesign <- data.frame(
  sampleID = colnames(quant_data$abundance)[-1],
  condition = c("Body","Body","Head","Head"))

#------------------------------------------------------------------
# Deal with the gff file 

# BASH:
# grep "exon" Dcitr_OGSv3.0_beta.gff3 > exons
# sed 's/ID=exon:/gene_id "/g' exons > exons1
# sed 's/\.1.*;Parent=/\.1"; transcript_id "/1' exons1 > exons2
# sed 's/;Name=.*$/"/g' exons2 > Dcitr_OGSv3_beta_dexseq.gtf

# Need to fix the double transcripts
# grep -v "," Dcitr_OGSv3_beta_dexseq.gtf > no_doubles
# grep "," Dcitr_OGSv3_beta_dexseq.gtf > doubles
# sed 's/,.*"/"/1' doubles > doubles_first_one_only

# R:
#doubles <- read_delim("doubles", "\t", escape_double = FALSE, 
#                      col_names = FALSE, trim_ws = TRUE)
#library(tidyverse)
#doubles_split <- doubles %>% 
#  mutate(X9=strsplit(X9, ",")) %>% 
#  unnest(X9)
#write.table(doubles_split, file="double_split", quote=F, col.names = F, row.names = F, sep = "\t")

# BASH:
# grep "gene_id" double_split > done
# sed 's/$/"/g' done > done1
# cat no_doubles done1 > working_file
# grep -v "gene_id" double_split > double_split1
# sed 's/Dcitr/transcript_id "Dcitr/g' double_split1 > double_split2

# sed 's/Dcitr.*"//g' double_split1 > double_split2
# awk -F '\t' 'BEGIN{OFS = "\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9$9}' double_split2 > double_split3
# sed 's/transcript_id/gene_id/1' double_split3 > double_split4
# sed 's/\.1\..*transcript_id/\.1"; transcript_id/g' double_split4 > double_split5
# cat working_file double_split5 > Dcitr_OGSv3_beta_dexseq2.gtf

# Import annotation
aSwitchList <- importRdata(
  isoformCountMatrix   = quant_data$counts,
  isoformRepExpression = quant_data$abundance,
  designMatrix         = myDesign,
  isoformExonAnnoation = "/Users/holliemarshall/Dropbox/Leicester_postdoc/Projects/lice/JCVI_LOUSE_1.0_genomic.gtf",
  showProgress = TRUE)

# Filter out genes with low expression and with no isoforms
SwitchListFiltered <- preFilter(
  switchAnalyzeRlist = aSwitchList,
  geneExpressionCutoff = 0,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE,
  dIFcutoff = 0.25) 

# Stopped here as only 4 isoforms being tested... what has gone wrong???







# FILTERED:
# 82 filtered for low expression and 17030 filtered because no annotated isoforms!
# 1419 left to test

#------------------------------------------------------------------
# Analyse remaining isoforms with DEXSeq for differntial usage
SwitchListAnalyzed <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = SwitchListFiltered,
  reduceToSwitchingGenes=TRUE)

extractSwitchSummary(SwitchListAnalyzed)
# OUTPUT:
#Comparison nrIsoforms nrSwitches nrGenes
#female vs male         37         26      21

switchingIso <- extractTopSwitches( 
  SwitchListAnalyzed, 
  filterForConsequences = F, 
  n = NA,              # n=NA: all features are returned
  extractGenes = T,    # when FALSE isoforms are returned
  sortByQvals = TRUE)

write.table(file="list_diff_alt_spliced_genes.txt",switchingIso$gene_id,
            sep = "\t", quote = F, col.names = T, row.names = F)

#------------------------------------------------------------------
# Make a plot for all of them as only 21 anyway
SwitchListAnalyzed[["isoformFeatures"]][["gene_name"]] <- SwitchListAnalyzed[["isoformFeatures"]][["gene_id"]]

# Can't get it to plot for some reason, maybe we just leave out the alternative splicing anyway as there are so few genes
switchPlotTopSwitches(
  switchAnalyzeRlist = SwitchListAnalyzed, 
  n = 1,
  filterForConsequences = FALSE, 
  splitFunctionalConsequences = F)



