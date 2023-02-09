# ------------------------------------------------------------
# Using parental RNA-Seq to confirm alt heterozygous snps
# ------------------------------------------------------------

setwd("~/Dropbox/Leicester_postdoc/Projects/lice/snp_counts_parents")

library(readr)
library(data.table)
library(dplyr)
library(doBy)
library(reshape2)

# Read in all data
file.list = list.files(("./"),pattern="*_snp_counts.txt")

read_file1 <- function(x){
  read_delim(x, "\t", escape_double = F, col_names = F, trim_ws = T)
}

samples <- lapply(file.list, read_file1)
sample_names <- list("BB_1","BB_2","HH_1","HH_2")
names(samples) <- sample_names

# Edit dataframes and filter by coverage
for(i in seq_along(samples)){
  samples[[i]] <- samples[[i]][,c(1,2,3,4,6,7,8,9)]
  colnames(samples[[i]]) <- c("chr", "SNP", "base", "coverage","A","C","G","T")
  samples[[i]] <- samples[[i]][samples[[i]]$coverage > 5,]
  samples[[i]]$A <- substring(samples[[i]]$A, 3)
  samples[[i]]$C <- substring(samples[[i]]$C, 3)
  samples[[i]]$G <- substring(samples[[i]]$G, 3)
  samples[[i]]$T <- substring(samples[[i]]$T, 3)
  samples[[i]]$A <- gsub(":.*", "",samples[[i]]$A)
  samples[[i]]$C <- gsub(":.*", "",samples[[i]]$C)
  samples[[i]]$G <- gsub(":.*", "",samples[[i]]$G)
  samples[[i]]$T <- gsub(":.*", "",samples[[i]]$T)
  samples[[i]]$base <- toupper(samples[[i]]$base)
  samples[[i]]$A <- as.numeric(samples[[i]]$A)
  samples[[i]]$C <- as.numeric(samples[[i]]$C)
  samples[[i]]$G <- as.numeric(samples[[i]]$G)
  samples[[i]]$T <- as.numeric(samples[[i]]$T)
}

# sum the info for both replicates
bb <- bind_rows(samples[[1]], samples[[2]]) %>%
  group_by(chr,SNP,base) %>%
  summarise_all(funs(sum(., na.rm = TRUE)))

hh <- bind_rows(samples[[3]], samples[[4]]) %>%
  group_by(chr,SNP,base) %>%
  summarise_all(funs(sum(., na.rm = TRUE)))

both <- list(bb, hh)

# Keep only SNPs where >90% reads come from that SNP
for(i in seq_along(both)){
  both[[i]]$prop_A <- both[[i]]$A/both[[i]]$coverage
  both[[i]]$prop_C <- both[[i]]$C/both[[i]]$coverage
  both[[i]]$prop_G <- both[[i]]$G/both[[i]]$coverage
  both[[i]]$prop_T <- both[[i]]$T/both[[i]]$coverage
  both[[i]] <- both[[i]][,c(1,2,9:12)]
  both[[i]] <- melt(both[[i]], id.vars=c("chr","SNP"))
  both[[i]] <- both[[i]][both[[i]]$value>0.9,]
  both[[i]]$variable <- gsub("prop_", "",both[[i]]$variable)
  both[[i]] <- both[[i]][,-4]
}


colnames(both[[1]])[3] <- "bb"
colnames(both[[2]])[3] <- "hh"
bb <- both[[1]] # 55185
hh <- both[[2]] # 59489

all_snps <- merge(bb, hh, by =c("chr","SNP")) # 50331
confirmed_diff_snps <- all_snps[!(all_snps$bb == all_snps$hh),] # 134

write.table(confirmed_diff_snps, file="snps_confirmed_with_RNA.txt", sep="\t", quote = F,
            row.names = F, col.names = T)
