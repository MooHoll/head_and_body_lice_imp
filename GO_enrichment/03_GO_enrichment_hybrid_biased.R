#-----------------------------------------------
# Script created by Alun Jones, see paper Bebane et al. (2019) Neonics and bumblebees...
#-----------------------------------------------

setwd("~/Dropbox/Leicester_postdoc/Projects/lice/GO_analysis/imprinted_analysis")
library(GOstats)
library(GSEABase)
library(treemap)

#-----------------------------------------------
# Read in background GO set and make compatible with GOstats

GO_annotations <- read.table("./background_all_hybrid_rnaseq_genes.txt")

#-----------------------------------------------
GO_annotations[,3] <- paste("IEA")
names(GO_annotations) <- c("GOIds","genes","evi")
GO_annotations[,3] <- paste("IEA")
GO_annotations <- GO_annotations[c(1,3,2)]
GO_annotations <- GO_annotations[-1,] # remove the original colnames

#-----------------------------------------------
# Create necessary objects

GO_frame <- GOFrame(GO_annotations,organism = "P.humanus")
goAllFrame <- GOAllFrame(GO_frame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
universe <- as.vector(unique(GO_annotations[,3]))

#-----------------------------------------------
# Read in gene's of interest 

#my_genes <- read.table("./body_expressed_all.txt", header = T)
#my_genes <- read.table("./head_expressed_all.txt", header = T)
#my_genes <- read.table("./maternally_expressed_all.txt", header = T)

#my_genes <- read.table("./body_expressed_limited.txt", header = T)
#my_genes <- read.table("./head_expressed_limited.txt", header = T)
my_genes <- read.table("./maternally_expressed_limited.txt", header = T)

#-----------------------------------------------
my_genes <- as.data.frame(na.omit(my_genes$gene_id))
colnames(my_genes) <- "genes"
my_genes <- as.vector(my_genes[,1])

#-----------------------------------------------
# Keep only genes with annotated GOs

my_genes <- my_genes[my_genes %in% universe]
length(my_genes)

# 6/26 body lice all
# 10/29 head lice all
# 6/18 maternal lice all

# 4/19 body lice limited
# 5/17 head lice limited
# 2/6 maternal lice limited

#-----------------------------------------------
# Set up paramters for hypergeometric test

Get_GO_params_all <- function(genes_of_i,universe,pvalue_cut){
  onto_terms <- c("BP","CC","MF")
  directions <- c("over","under")
  param_list <- list()
  name_1 <- list()
  for(i in 1:3){
    for(j in 1:2){
      name_1 <- c(name_1,paste(onto_terms[i],directions[j],sep = "_"))
      parameters <- GSEAGOHyperGParams(name="Hygeo params",
                                       geneSetCollection = gsc,
                                       universeGeneIds = universe,
                                       geneIds = my_genes,
                                       ontology = paste(onto_terms[i]),
                                       pvalueCutoff = pvalue_cut,
                                       conditional = T,testDirection = paste(directions[j]))
      param_list <- c(param_list,parameters)
    }
  }
  names(param_list) <- name_1
  return(param_list)
}

param_list <- Get_GO_params_all(genes_of_i = DE_Genes_A,universe = universe,
                                pvalue_cut = 0.05)


#-----------------------------------------------
# Hypergeometric test

Hyper_G_test <- function(param_list){
  Hyper_G_list <- list()
  for(i in 1:length(param_list)){
    res <- hyperGTest(param_list[[i]])
    Hyper_G_list <- c(Hyper_G_list,res)
  }
  names(Hyper_G_list) <- names(param_list)
  return(Hyper_G_list)
}


GO_enrichment <- Hyper_G_test(param_list = param_list)


#-----------------------------------------------
# Choose what you want to look for, here the choice is biological process over-represented 

Result <- summary(GO_enrichment[["BP_over"]])

#-----------------------------------------------
# FDR correction

Result_FDR <- Result[p.adjust(Result$Pvalue,method = "fdr") < 0.05,]

#-----------------------------------------------
# Make an output for REVIGO and write out

REVIGO <- Result_FDR[,1:2]

#write.table(REVIGO,"./outputs/body_expressed_all_enriched_GOs.txt",row.names = F,sep = "\t",quote = F)
#write.table(REVIGO,"./outputs/head_expressed_all_enriched_GOs.txt",row.names = F,sep = "\t",quote = F)
#write.table(REVIGO,"./outputs/maternal_expressed_all_enriched_GOs.txt",row.names = F,sep = "\t",quote = F)

#write.table(REVIGO,"./outputs/body_expressed_limited_enriched_GOs.txt",row.names = F,sep = "\t",quote = F)
#write.table(REVIGO,"./outputs/head_expressed_limited_enriched_GOs.txt",row.names = F,sep = "\t",quote = F)
write.table(REVIGO,"./outputs/maternal_expressed_limited_enriched_GOs.txt",row.names = F,sep = "\t",quote = F)

