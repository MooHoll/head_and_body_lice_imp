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

# Import annotation
aSwitchList <- importRdata(
  isoformCountMatrix   = quant_data$counts,
  isoformRepExpression = quant_data$abundance,
  designMatrix         = myDesign,
  isoformExonAnnoation = "/Users/holliemarshall/Dropbox/Leicester_postdoc/Projects/lice/genome/JCVI_LOUSE_1.0_genomic.gtf",
  showProgress = TRUE)

# Filter out genes with low expression and with no isoforms
SwitchListFiltered <- preFilter(
  switchAnalyzeRlist = aSwitchList,
  geneExpressionCutoff = 0,
  isoformExpressionCutoff = 0,
  removeSingleIsoformGenes = TRUE,
  dIFcutoff = 0.10) # Filtering only leaves 4 isoforms! Why!


# Analyse remaining isoforms with DEXSeq for differntial usage
SwitchListAnalyzed <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = aSwitchList,
  reduceToSwitchingGenes=TRUE) # This gives us no genes at all significant







extractSwitchSummary(SwitchListAnalyzed)

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



