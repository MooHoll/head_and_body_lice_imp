###-------------------------------------------------------------------
# Run imprinting model
###-------------------------------------------------------------------

setwd("~/Dropbox/Leicester_postdoc/Projects/lice/snp_counts")

library(readr)
library(car)
library(lsmeans)
library(reshape2)
library(ggplot2)

###-------------------------------------------------------------------
# Read in data
all_data <- read_delim("all_data_for_imp_model.txt", 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE)
all_data$gene_id <- as.factor(all_data$gene_id)
all_data$cross <- as.factor(all_data$cross)

###-------------------------------------------------------------------
# Test with one gene
one_gene <- all_data[all_data$gene_id=="gene-Phum_PHUM000210",]

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
model_output$threshold[model_output$avgpropmatexpr>0.6] <- "maternal"
model_output$threshold[model_output$avgpropmatexpr<0.4] <- "paternal"
model_output$actually_sig <- "no"
model_output$actually_sig[model_output$sig_in_both == 1 & !model_output$threshold == "unbiased"] <- "yes"

# Make a dataframe for plotting
for_plot <- model_output[,c(1,2,9,15)]
for_plot1 <- dcast(for_plot, gene_id + actually_sig ~ direction_cross, value.var = "avgpropmatexpr")

# Make a plot
ggplot(for_plot1, aes(x=hb, y=bh, colour =actually_sig))+
  geom_point()

