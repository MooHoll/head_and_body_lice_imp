#------------------------------------------------------------------
# Making a file with FPKM and log2FC for later incorporation with methylation data
#------------------------------------------------------------------

setwd("~/Dropbox/Leicester_postdoc/Projects/lice/RNA_seq_pure/counts/Pure_data")

library(readr)
library(reshape2)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggforce)

#------------------------------------------------------------------
# log2FC from diff exp
diff_exp_log2FC_all <- read_delim("diff_exp_output_log2FC_all_genes.txt", 
                                               "\t", escape_double = FALSE, trim_ws = TRUE)
diff_exp_log2FC_all <- diff_exp_log2FC_all[,c(2,5,6,7)]
colnames(diff_exp_log2FC_all)[4] <- "gene_id"

#------------------------------------------------------------------
# all samples for FPKM
BB1 <- read_delim("BB_M_1.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
BB1 <- BB1[,c(1,7)]
colnames(BB1) <- c("gene_id", "BB1_FPKM")

BB2 <- read_delim("BB_M_2.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
BB2 <- BB2[,c(1,7)]
colnames(BB2) <- c("gene_id", "BB2_FPKM")

HH1 <- read_delim("HH_M_1.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
HH1 <- HH1[,c(1,7)]
colnames(HH1) <- c("gene_id", "HH1_FPKM")

HH2 <- read_delim("HH_M_2.genes.results", "\t", escape_double = FALSE, trim_ws = TRUE)
HH2 <- HH2[,c(1,7)]
colnames(HH2) <- c("gene_id", "HH2_FPKM")

#------------------------------------------------------------------
# Put it all together
fpkm_data <- Reduce(merge, list(BB1,BB2,HH1,HH2))

fpkm_data$BB_fpkm_mean <- (fpkm_data$BB1_FPKM + fpkm_data$BB2_FPKM)/2
fpkm_data$HH_fpkm_mean <- (fpkm_data$HH1_FPKM + fpkm_data$HH2_FPKM)/2

all_data <- merge(fpkm_data, diff_exp_log2FC_all, by = "gene_id") #9830

# Measure of specificity for gene expression
all_data$BB_sq <- all_data$BB_fpkm_mean*all_data$BB_fpkm_mean
all_data$HH_sq <- all_data$HH_fpkm_mean*all_data$HH_fpkm_mean
all_data$SPM <- all_data$BB_sq/(all_data$BB_sq +all_data$HH_sq )

ggplot(all_data, aes ( x= SPM))+
  geom_histogram(colour="black", bins=50)+
  xlab("SPM Relative to Body Lice")+
  ylab("Number of Genes")+
  theme_bw()+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        legend.text = element_text(size=22))

#------------------------------------------------------------------
# Add gene categories for plotting
all_data$diff_exp <- "no"
all_data$diff_exp[all_data$padj < 0.05] <- "yes"
all_data$diff_exp[all_data$log2FoldChange < 1.5 & all_data$log2FoldChange > -1.5] <- "no"

table(all_data$diff_exp) # 47 = yes

all_data$category <- "unbiased"
all_data$category[all_data$diff_exp =="yes" & all_data$log2FoldChange > 1.5] <- "HH_biased" # 13
all_data$category[all_data$diff_exp =="yes" & all_data$log2FoldChange < -1.5] <- "BB_biased" # 33
all_data$category[all_data$diff_exp =="yes" & all_data$log2FoldChange > 10] <- "HH_biased_extreme" # 0
all_data$category[all_data$diff_exp =="yes" & all_data$log2FoldChange < -10] <- "BB_biased_extreme" # 1

table(all_data$category)

all_data$category[all_data$diff_exp =="yes" & all_data$BB_fpkm_mean==0] <- "HH_limited" # 0
all_data$category[all_data$diff_exp =="yes" & all_data$HH_fpkm_mean==0] <- "BB_limited" # 2

table(all_data$category)

# Take an eyeball of the male limited genes
head(all_data)
look <- all_data[all_data$category=="BB_limited",]
boxplot(look$BB_fpkm_mean)
mean(look$BB_fpkm_mean)
median(look$BB_fpkm_mean)
range(look$BB_fpkm_mean)
nrow(look[look$BB_fpkm_mean<10,]) #Hmmm 1/2, not too sure about this

# Goodness of fit
observed = c(34, 13)    # observed frequencies 
expected = c(0.5, 0.5)    # expected proportions

chisq.test(x = observed,
           p = expected) 
#X-squared = 9.383, df = 1, p-value = 0.00219 - sig more BB overexpressed compared to HH


write.table(all_data, file="all_gene_expression_data.txt", sep="\t", quote=F, col.names = T, row.names = F )
#------------------------------------------------------------------
# Plots
all_data_plot <- all_data[!all_data$category == "unbiased",]
all_data_plot$biased[all_data_plot$diff_exp =="yes" & all_data_plot$log2FoldChange > 1.5] <- "HH_biased"
all_data_plot$biased[all_data_plot$diff_exp =="yes" & all_data_plot$log2FoldChange < -1.5] <- "BB_biased"
all_data_plot$category <- gsub(".*_","",all_data_plot$category)

ggplot(all_data_plot, aes(x=biased, fill=category))+
  geom_bar()+
  theme_bw()+
  xlab("Gene Expression Category")+
  ylab("Number of Genes")+
  scale_fill_manual("",breaks=c("biased", "extreme","limited"),
                    labels = c("Fold-change > 1.5","Fold-change > 10","Ecotype Limited"),
                    values = c("grey75", "grey45","grey14"))+
  theme_bw()+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        legend.text = element_text(size=22))
#------------------------------------------------------------------
head(all_data)
all_data$BB_log_fpkm <- log10(all_data$BB_fpkm_mean)
all_data$HH_log_fpkm <- log10(all_data$HH_fpkm_mean)

scatter_data <- all_data[,c(15:17)]

ggplot(scatter_data, aes(x=HH_log_fpkm, y=BB_log_fpkm, colour=category))+
  geom_point(size=3)+
  xlim(0,4)+
  ylim(0,4)+
  scale_colour_manual("",breaks=c("unbiased","BB_biased","BB_limited","HH_biased"),
                      labels=c("Unbiased","BB Biased","BB Limited","HH Biased"),
                      values=c("grey","#009E73","darkgreen","#0072B3"))+
  theme_bw()+
  xlab("Log10(HH FPKM)")+
  ylab("Log10(BB FPKM)")+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        legend.text = element_text(size=22))
