# ------------------------------------------------------------
# Using parental RNA-Seq to confirm alt homozygous snps and annotate with gene name
# ------------------------------------------------------------

setwd("~/Dropbox/Research/Leicester_postdoc/Projects/lice/snp_counts")

library(readr)
library(data.table)
library(sqldf)
library(doBy)
library(foreach)
library(doParallel)
library(reshape2)

# Read in all data
file.list = list.files(("./"),pattern="*counts.txt")

read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = F, trim_ws = T)
}

samples <- lapply(file.list, read_file1)
sample_names <- list("bh_1","bh_4","bh_5","bh_6",
                     "hb_1","hb_2","hb_3","hb_5","hb_6")
names(samples) <- sample_names

# Edit dataframes, coverage > 10 per SNP and make final columns useable
for(i in seq_along(samples)){
  samples[[i]] <- samples[[i]][,c(1,2,3,4,6,7,8,9)]
  colnames(samples[[i]]) <- c("chr", "SNP", "base", "coverage","A","C","G","T")
  samples[[i]] <- samples[[i]][samples[[i]]$coverage > 10,]
  samples[[i]]$A <- substring(samples[[i]]$A, 3)
  samples[[i]]$C <- substring(samples[[i]]$C, 3)
  samples[[i]]$G <- substring(samples[[i]]$G, 3)
  samples[[i]]$T <- substring(samples[[i]]$T, 3)
  samples[[i]]$A <- gsub(":.*", "",samples[[i]]$A)
  samples[[i]]$C <- gsub(":.*", "",samples[[i]]$C)
  samples[[i]]$G <- gsub(":.*", "",samples[[i]]$G)
  samples[[i]]$T <- gsub(":.*", "",samples[[i]]$T)
  samples[[i]]$base <- toupper(samples[[i]]$base)
  samples[[i]] <- samples[[i]][,-c(3,4)]
  samples[[i]]$sample <- paste(names(samples[i]))
}

# Read in actual called snps
bb_snps <- read_delim("../bl_snps.txt", delim = "\t", 
                      escape_double = FALSE, col_names = FALSE, 
                      trim_ws = TRUE)
colnames(bb_snps)<-c("chr","SNP","hh","bb") 

hh_snps <- read_delim("../hl_snps.txt", delim = "\t", 
                      escape_double = FALSE, col_names = FALSE, 
                      trim_ws = TRUE)
colnames(hh_snps)<-c("chr","SNP","bb","hh")

# Check for any positions where there is a SNP in both the head and body louse at the same position
look <- merge(bb_snps, hh_snps, by = c("chr","SNP")) # 0! ok happy days

all_snps <- rbind(hh_snps, bb_snps)


# Just use the confirmed SNPs first
#snps_confirmed_with_RNA <- read_delim("~/Dropbox/Leicester_postdoc/Projects/lice/snp_counts_parents/snps_confirmed_with_RNA.txt", 
#                                        delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# Merge with SNPs that are confirmed and rename columns by matching SNPs
for(i in seq_along(samples)){
  samples[[i]] <- merge(samples[[i]], all_snps, by=c("chr","SNP"))
  samples[[i]]$bb_count <- samples[[i]][-ncol(samples[[i]])][cbind(seq_len(nrow(samples[[i]])),
                                           match(samples[[i]]$bb, names(samples[[i]])[-ncol(samples[[i]])]))]
  samples[[i]]$hh_count <- samples[[i]][-ncol(samples[[i]])][cbind(seq_len(nrow(samples[[i]])),
                                           match(samples[[i]]$hh, names(samples[[i]])[-ncol(samples[[i]])]))]
  samples[[i]] <- samples[[i]][,c(1,2,7,10,11)]
  samples[[i]]$cross <- paste(substr(samples[[i]]$sample, 1, 2))
}

# ------------------------------------------------------------------
# Annotate with gene
louse_genes <- read_delim("~/Dropbox/Leicester_postdoc/Projects/lice/genome/louse_genes.bed", 
                          delim = "\t", escape_double = FALSE, 
                          col_names = FALSE, trim_ws = TRUE)
colnames(louse_genes) <- c("chr","start","end","gene_id")

registerDoParallel(cores = 2)
foreach(i = seq_along(samples)) %dopar% {
  df <- samples[[i]]
  output <- sqldf("SELECT sample.chr,
                    sample.SNP,
                    sample.sample,
                    sample.bb_count,
                    sample.hh_count,
                    sample.cross,
                    genes.chr,
                    genes.start,
                    genes.end,
                    genes.gene_id
                    FROM df AS sample
                    LEFT JOIN louse_genes AS genes
                    ON sample.chr = genes.chr
                    AND (sample.SNP >= genes.start AND sample.SNP <= genes.end)")
  output <- output[!is.na(output$gene_id),]
  output <- output[,c(3,4,5,6,10)]
  output <- melt(output, id.vars=c("sample","cross","gene_id"))
  colnames(output) <- c("sample","cross","gene_id","line_origin","count")
  output$count <- as.numeric(output$count)
  check <- summaryBy(count ~ sample + cross + gene_id + line_origin, data=output, FUN=sum) 
  myfile <- file.path("./", paste0(names(samples[i]),"_","counts_by_genes.txt"))
  write.table(check, file=myfile, quote=F, sep="\t", row.names=F)
}

