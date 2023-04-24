#---------------------------------------------------------
# Expression levels of DNMT genes
#---------------------------------------------------------

setwd("~/Dropbox/Leicester_postdoc/Projects/lice/RNA_seq_pure/counts/Pure_data")

library(readr)
library(reshape2)
library(ggplot2)
#---------------------------------------------------------

all_gene_expression_data <- read_delim("all_gene_expression_data.txt", 
                                       "\t", escape_double = FALSE, trim_ws = TRUE)

# DNMT1 gene-Phum_PHUM556380
# DNMT3a gene-Phum_PHUM041250

head(all_gene_expression_data)

dnmt1_1 <- all_gene_expression_data[all_gene_expression_data$gene_id=="gene-Phum_PHUM556380",]
# Unbiased and low exp

dnmt3a <- all_gene_expression_data[all_gene_expression_data$gene_id=="gene-Phum_PHUM041250",]
# Unbiased and slightly higher exp

dnmt1_1 <- dnmt1_1[,c(1:5)]
dnmt1_1 <- melt(dnmt1_1)
dnmt1_1$Ecotype <- c("BB","BB","HH","HH")

ggplot(dnmt1_1, aes(x=Ecotype, y=value, fill=Ecotype))+
  geom_boxplot()+
  geom_jitter(size=5)+
  theme_bw()+
  xlab("Ecotype")+
  ylab("Expression Level (FPKM)")+
  ggtitle("DNMT1")+
  scale_fill_manual("", breaks=c("BB","HH"),
                      values = c("#009E73","#0072B3"))+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        legend.position = "none",
        title = element_text(size=22))


dnmt3a <- dnmt3a[,c(1:5)]
dnmt3a <- melt(dnmt3a)
dnmt3a$Ecotype <- c("BB","BB","HH","HH")

ggplot(dnmt3a, aes(x=Ecotype, y=value, fill=Ecotype))+
  geom_boxplot()+
  geom_jitter(size=5)+
  theme_bw()+
  xlab("Ecotype")+
  ylab("Expression Level (FPKM)")+
  ggtitle("DNMT3a")+
  scale_fill_manual("", breaks=c("BB","HH"),
                    values = c("#009E73","#0072B3"))+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        legend.position = "none",
        title = element_text(size=22))
