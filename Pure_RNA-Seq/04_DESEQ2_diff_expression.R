#------------------------------------------------------------------
# Differential Gene Expression
#------------------------------------------------------------------

setwd("~/Dropbox/Leicester_postdoc/Projects/lice/RNA_seq_pure/counts/Pure_data")
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
library(PoiClaClu)
library(apeglm)
library(tximportData)
library(tximport)

#------------------------------------------------------------------
# Make sample metadata
sample <- c("BB_1","BB_2", "HH_1","HH_2")
ecotype <- c("Body","Body","Head","Head")
file_name <- list.files("./", pattern="*.genes.results")
sample_info <- data.frame(sample, ecotype, file_name)
rownames(sample_info) <- sample_info$sample

# Read in files using tximport
files <- file.path("./", paste0(sample_info$file_name))
names(files) <- paste0(sample_info$sample)
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
head(txi.rsem$counts)
txi.rsem$length[txi.rsem$length == 0] <- 1

# Make deseq2 object
dds <- DESeqDataSetFromTximport(txi.rsem,
                                   colData = sample_info,
                                   design = ~ ecotype)

#------------------------------------------------------------------
# remove features with low counts
dds = dds[ rowMeans(counts(dds)) > 10, ] 
nrow(dds) # 9830/10992

# rlog transform counts (by average Tx length and correcting for library size)
rld = rlog(dds, blind=FALSE)

#------------------------------------------------------------------
# PCA plot
data = plotPCA(rld, intgroup = c("ecotype"), returnData=TRUE)
percentVar = round(100 * attr(data, "percentVar"))

ggplot(data, aes(PC1, PC2, color=ecotype)) + geom_point(size=14) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #coord_fixed(ratio=4,clip = "on")+
  geom_text_repel(aes(label=name), size=10,show.legend=FALSE, 
                  point.padding = 2, box.padding = 1,
                  segment.color = 'transparent') +
  scale_colour_manual("", breaks=c("Body","Head"),
                      values = c("#009E73","#0072B3"))+
  theme_bw()+
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=24),
        legend.text=element_text(size=24))

# two 1st samples plotted against each other to check consistency (for rlog and log2) 
par( mfrow = c( 1, 2 ) )
dds = estimateSizeFactors(dds)
plot(log2(counts(dds, normalized=TRUE)[,c(1,2)] + 1),
     pch=16, cex=0.3, cex.lab=1.5, cex.axis=1.5)
plot(assay(rld)[,c(1,2)],
     pch=16, cex=0.3, cex.lab=1.5, cex.axis=1.5)
par( mfrow = c( 1, 1 ) )

# estimate size factors = normalize for dispersion
dds = DESeq2::estimateSizeFactors(dds)
dds = estimateDispersions(dds)
plotDispEsts(dds, xlab= "Mean of Normalised Counts",
             ylab= "Dispersion", cex=1.0, cex.lab=1.45, cex.axis=1.45)

# check sample distances
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$dex, rld$cell, sep = " - " )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         show_colnames = T,
         fontsize = 18)

# check sample distances using the poisson method
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( rld$dex, rld$cell, sep=" - " )
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors,
         show_colnames = T,
         fontsize = 18)

#------------------------------------------------------------------
# differential expression (used Benjamini-Hochberg adjustment and a p-val of <0.1)
dds = DESeq2::DESeq(dds, parallel=TRUE)
resultsNames(dds)

# Head = +Log2FC and Body = -Log2FC
res=results(dds, name="ecotype_Head_vs_Body")
summary(res)

#distribution of coefficents of the model
plotMA(res, ylim=c(-5,5),cex=1.0, cex.lab=1.45, 
       cex.axis=1.45, xlab='Mean of Normalised Counts',
       ylab='Log2 Fold Change')

# plot of p-vals excluding genes with very small counts
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white", xlab= "P-value", ylab="Frequency",
     cex.axis=1.45, cex.lab=1.45, main="")

# Order the results by fold change to see the biggest changing genes
res_ordered=res[order(res$log2FoldChange),]
head(res_ordered)
res_ordered<-as.data.frame(res_ordered)
res_ordered$gene<-rownames(res_ordered)
rownames(res_ordered)<-c()
write.table(as.data.frame(res_ordered), file="diff_exp_output_log2FC_all_genes.txt",
            sep="\t", quote = F, col.names = T, row.names = F)

#out of 9830 with nonzero total read count
res_significant<- subset(res, log2FoldChange > 1.5 | log2FoldChange < -1.5)
res_significant <- subset(res_significant, padj < 0.05)
nrow(res_significant) #47 total differentially expressed genes
nrow(res_significant[res_significant$log2FoldChange > 1.5,]) #13 upregulated in head
nrow(res_significant[res_significant$log2FoldChange < -1.5,]) #34 upregulated in body

#------------------------------------------------------------------
# heatmap of top differentially expressed
n=47
topdiff = head(c(1:nrow(res))[order(res$padj)],n)

my_colors = list(
  ecotype = c(Body = "#009E73", Head ="#0072B3"))


mat = assay(rld)[ topdiff, ]
mat = mat - rowMeans(mat)
df2 = as.data.frame(colData(rld)[,c("ecotype"),drop=FALSE])
colnames(df2)<-c("ecotype")

pheatmap(mat, annotation_col=df2,
         show_rownames = F,
         fontsize = 16,
         annotation_colors = my_colors)

#------------------------------------------------------------------
# Plot diff expressed
n=12
selGenes = head(rownames(res)[order(res$padj)],n)
data = do.call(rbind, lapply(selGenes, function(gene) data.frame(gene=gene, 
                                                                 plotCounts(dds, gene=gene, intgroup=c("ecotype"), returnData=TRUE))))
ggplot(data, aes(x=ecotype, y=count, fill=ecotype)) + 
  geom_boxplot(outlier.color="black", position=position_dodge(width=0.7), width=0.5) +
  facet_wrap(~gene) +
  xlab("Ecotype") + ylab("Normalised read count") + 
  scale_y_log10() + 
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.title = element_text(size=20),
        legend.text=element_text(size=20),
        axis.text.y=element_text(size=18),
        plot.title = element_text(size=22))+
  scale_fill_manual("", breaks=c("Body","Head"),
                      values = c("#009E73","#0072B3"))

#------------------------------------------------------------------
# Volcano plot
res_df <- as.data.frame(res)
res_df$gene <- row.names(res_df)

res_significant_df <- as.data.frame(res_significant)
res_significant_df$gene <- row.names(res_significant_df)

res_df$sig <- "no"
res_df$sig[res_df$gene %in% res_significant_df$gene] <- "yes"

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), colour = sig)) + 
  geom_point(size=1.5)+
  scale_colour_manual("", breaks=c("no","yes"),
                      values = c("black","red"))+
  xlim(-10,10)+
  theme_bw()+
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.position = "none")

#------------------------------------------------------------------
# Pull out a list of just the differentially expressed genes for the supplementary



