###-------------------------------------------------------------------
# Run imprinting model
###-------------------------------------------------------------------

setwd("~/Dropbox/Leicester_postdoc/Projects/lice/snp_counts")

library(readr)
library(car)
library(lsmeans)
library(reshape2)
library(ggplot2)
library(tidyr)


###-------------------------------------------------------------------
# Read in data
all_data <- read_delim("all_data_for_imp_model.txt", 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE)
all_data$gene_id <- as.factor(all_data$gene_id)
all_data$cross <- as.factor(all_data$cross)

###-------------------------------------------------------------------
# Test with one gene
one_gene <- all_data[all_data$gene_id=="gene-Phum_PHUM006440",]

fit <- glm( cbind(maternal, paternal) ~ cross, 
           family = quasibinomial(link=logit), data=one_gene)
summary(fit) # Wald z tests
Anova(fit,type="III")

model_output <- data.frame(summary(suppressMessages(lsmeans(fit,~ cross))))
model_output$z <- model_output$lsmean/model_output$SE # z scores for deviation from 50:50
model_output$p <- 2*pnorm(abs(model_output$z),lower.tail = FALSE) # 2 sided p values for deviation from 50:50 (not corrected for mult testing)
model_output$lsmean <- plogis(model_output$lsmean) # backtransform to original scale
model_output$asymp.LCL <- plogis(model_output$asymp.LCL)
model_output$asymp.UCL <- plogis(model_output$asymp.UCL)


###-------------------------------------------------------------------
# Model for all genes
model_output = do.call(rbind,
                           lapply(levels(all_data$gene_id),
                                    function (gene) { 
                                      dfsubs = all_data[all_data$gene_id==gene,]
                                      
                                      fit = glm( cbind(maternal, paternal) ~ cross, 
                                                 family = quasibinomial(link=logit), data=dfsubs)

                                      model_output = data.frame(summary(suppressMessages(lsmeans(fit,~ cross))))
                                      model_output$z = model_output$lsmean/model_output$SE # z scores for deviation from 50:50
                                      model_output$p = 2*pnorm(abs(model_output$z),lower.tail = FALSE) # 2 sided p values for deviation from 50:50 (not corrected for mult testing)
                                      model_output$lsmean = plogis(model_output$lsmean) # backtransform to original scale
                                      model_output$asymp.LCL = plogis(model_output$asymp.LCL)
                                      model_output$asymp.UCL = plogis(model_output$asymp.UCL)
                                      df=data.frame(gene_id=gene,model_output)
                                      df=df[,-5] # remove the df column
                                      colnames(df)=c("gene_id","direction_cross","proportion_mat_exp","SE","propmatexpr.LCL","propmatexpr.UCL","z","p")
                                      df$avgpropmatexpr=(df$proportion_mat_exp+df$propmatexpr.LCL+df$propmatexpr.UCL)/3
                                      df
                                    }
                           ) )

# Benjamini-Hochberg correction for multiple testing
model_output$padj <- p.adjust(model_output$p, method="BH") 

# Add a column showing the significant rows
model_output$sig <- (model_output$padj<0.05)*1

# Add a column showing if the expression bias is maternal (1) or paternal (0)
model_output$mat_exp <- (model_output$avgpropmatexpr>0.5)*1

# Label genes present only in both crosses
new <- model_output[,c(1,2,11)]
new1 <- dcast(new, gene_id ~ direction_cross)
new1$sig_in_both[new1$bh == 1 & new1$hb == 1] <- "yes"
model_output$sig_in_both <- (model_output$gene_id %in% new1$gene_id[new1$sig_in_both=="yes"])*1

# Hard filter on threshold
model_output$threshold <- "unbiased"
model_output$threshold[model_output$avgpropmatexpr>=0.6] <- "maternal"
model_output$threshold[model_output$avgpropmatexpr<=0.4] <- "paternal"
model_output$actually_sig <- "no"
model_output$actually_sig[model_output$sig_in_both == 1 & !model_output$threshold == "unbiased"] <- "yes"
table(model_output$actually_sig) # 147 genes... means at least one gene passed the sig threshold but is not within the hard threshold of 0.6/0.4

write.table(model_output, file="imp_model_output.txt", sep="\t", quote = F, col.names = T, row.names = F)


# be sure the hard threshold has worked
look <- model_output[model_output$threshold=="unbiased" & model_output$sig_in_both == 1,] # just one gene

# Add column to can colour the dots in the plot
for_plot <- model_output[,c(1,2,9,14,15)]
for_plot <- for_plot[for_plot$actually_sig =="yes",] #147 (two rows per gene = actually )
# throw out that naughty gene above
for_plot <- for_plot[!for_plot$gene_id %in% look$gene_id,] # 146, sorted
for_plot_wide <- dcast(for_plot, gene_id + actually_sig ~ threshold)
maternal_genes <- as.data.frame(for_plot_wide$gene_id[for_plot_wide$maternal==2]) # 18 genes
colnames(maternal_genes) <- "gene_id"
maternal_genes$lineage <- "maternal"

for_plot_lineage <- for_plot[!for_plot$gene_id %in% maternal_genes$gene_id,]

for_plot_lineage <- spread(for_plot_lineage, key = threshold, value = avgpropmatexpr )
for_plot_lineage$lineage <-"body"
for_plot_lineage$lineage[for_plot_lineage$direction_cross=="hb" &
                           !is.na(for_plot_lineage$maternal)] <-"head"
for_plot_lineage <- for_plot_lineage[!is.na(for_plot_lineage$maternal),]
for_plot_lineage <- for_plot_lineage[,c(1,6)]
table(for_plot_lineage$lineage) # 26 body biased, 29 head biased

# Put in the maternal data
sig_genes_with_where_from <- rbind(for_plot_lineage, maternal_genes) # 73 total

# Make a dataframe for plotting
for_plot <- model_output[,c(1,2,9,15)] # note: 482 genes total tested
for_plot <- dcast(for_plot, gene_id + actually_sig ~ direction_cross, value.var = "avgpropmatexpr")

for_plot <- merge(for_plot, sig_genes_with_where_from, all=T)
for_plot$lineage[is.na(for_plot$lineage)] <- "unbiased"
table(for_plot$lineage)
#body     head maternal unbiased 
#26       29       18      409 

# Make a plot
ggplot(for_plot, aes(x=hb, y=bh, colour =lineage))+
  geom_point(size=2.5)+
  xlim(0,1)+
  ylim(0,1)+
  theme_bw()+
  xlab("Proportion maternal expression (HB cross)")+
  ylab("Proportion maternal expression (BH cross)")+
  geom_hline(yintercept=0.4, linetype="dashed", color = "black")+
  geom_vline(xintercept=0.4, linetype="dashed", color = "black")+
  geom_hline(yintercept=0.6, linetype="dashed", color = "black")+
  geom_vline(xintercept=0.6, linetype="dashed", color = "black")+
  theme(axis.title = element_text(size=19),
        axis.text=element_text(size=18),
        legend.text = element_text(size=16))+
  scale_colour_manual(breaks=c("body","head", "maternal", "unbiased"),
                      labels=c("Body Louse Biased", "Head Louse Biased", "Materally Biased", "Unbiased"),
                    values = c("#009E73","#0072B3","deeppink","grey"),
                    name ="")

# How many are completely silenced?
head(for_plot)

silenced <- for_plot[for_plot$actually_sig=="yes" & (for_plot$bh>0.99 | for_plot$hb>0.99),]
for_plot$silenced <- "no"
for_plot$silenced[for_plot$gene_id %in% silenced$gene_id] <- "yes"
table(for_plot$lineage[for_plot$silenced=="yes" & for_plot$actually_sig=="yes"])
#  body     head maternal 
# 19       17        6

for_bar_plot <- for_plot[for_plot$gene_id %in% sig_genes_with_where_from$gene_id,]

ggplot(for_bar_plot, aes(x=lineage, fill=silenced))+
  geom_bar(stat="count")+
  theme_bw()+
  xlab("Expression Category")+
  ylab("Number of Genes")+
  scale_fill_manual("Expression Level",breaks=c("no", "yes"),
                    labels = c("Biased","Limited"),
                    values = c("grey45", "grey14"))+
  scale_x_discrete(breaks = c("body","head","maternal"), 
                   labels=c("Body\nLouse", "Head\nLouse", "Maternal"))+
  theme_bw()+
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=22),
        legend.text = element_text(size=22),
        legend.title = element_text(size=22))

# Write out a useable dataframe for overlaps and GO analysis etc
head(for_plot)

write.table(for_plot, file="~/Dropbox/Leicester_postdoc/Projects/lice/imp_genes_output.txt", sep="\t", quote = F, row.names = F, col.names = T)

# Gene lists for GO analysis
head(for_plot)

setwd("~/Dropbox/Leicester_postdoc/Projects/lice/GO_analysis/imprinted_analysis")

# all - both biased and limited
maternal_all <- for_plot[for_plot$actually_sig=="yes" & for_plot$lineage=="maternal",] #18
maternal_all <- as.data.frame(maternal_all$gene_id)
colnames(maternal_all) <- "gene_id"
write.table(maternal_all, file="maternally_expressed_all.txt", sep="\t", quote = F, col.names = T, row.names = F)

body_all <- for_plot[for_plot$actually_sig=="yes" & for_plot$lineage=="body",] #26
body_all <- as.data.frame(body_all$gene_id)
colnames(body_all) <- "gene_id"
write.table(body_all, file="body_expressed_all.txt", sep="\t", quote = F, col.names = T, row.names = F)

head_all <- for_plot[for_plot$actually_sig=="yes" & for_plot$lineage=="head",] #29
head_all <- as.data.frame(head_all$gene_id)
colnames(head_all) <- "gene_id"
write.table(head_all, file="head_expressed_all.txt", sep="\t", quote = F, col.names = T, row.names = F)

# and the limited ones
maternal_limited <- for_plot[for_plot$actually_sig=="yes" & for_plot$lineage=="maternal" & for_plot$silenced=="yes",] #6
maternal_limited <- as.data.frame(maternal_limited$gene_id)
colnames(maternal_limited) <- "gene_id"
write.table(maternal_limited, file="maternally_expressed_limited.txt", sep="\t", quote = F, col.names = T, row.names = F)

body_limited <- for_plot[for_plot$actually_sig=="yes" & for_plot$lineage=="body" & for_plot$silenced=="yes",] #19
body_limited <- as.data.frame(body_limited$gene_id)
colnames(body_limited) <- "gene_id"
write.table(body_limited, file="body_expressed_limited.txt", sep="\t", quote = F, col.names = T, row.names = F)

head_limited <- for_plot[for_plot$actually_sig=="yes" & for_plot$lineage=="head" & for_plot$silenced=="yes",] #17
head_limited <- as.data.frame(head_limited$gene_id)
colnames(head_limited) <- "gene_id"
write.table(head_limited, file="head_expressed_limited.txt", sep="\t", quote = F, col.names = T, row.names = F)

# Make final dataframe for supplementary
head(model_output)
head(for_plot)
for_supp <- for_plot[,c(1,5)]

both <- merge(model_output, for_supp) #481
write.table(both, file="imp_model_output_for_supp.txt", sep="\t", quote = F, col.names = T, row.names = F)
