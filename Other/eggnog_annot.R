#----------------------------------------------------------
# Make a Gene Ontology annotation file
#----------------------------------------------------------

# Made a protein fasta file from the genome on ALICE:
#conda activate agat_env
#agat_sp_extract_sequences.pl --gff JCVI_LOUSE_1.0_genomic.gff -f JCVI_LOUSE_1.0_genomic.fa -p -o louse_proteins.fa

# Uploaded the protein file to eggnog, standard parameters:
# http://eggnog-mapper.embl.de

# Downloaded the excel output
# Removed by hand all info from the .xlxs file except gene id and GO annotation

setwd("~/Dropbox/Leicester_postdoc/Projects/lice/genome")

library(tidyr)
library(dplyr)
library(readxl)

out_emapper_annotations <- read_excel("out.emapper.annotations.xlsx")
colnames(out_emapper_annotations) <- c("gene_id","GO_term")

new <- out_emapper_annotations %>% 
  mutate(GO_term = strsplit(as.character(GO_term), ",")) %>% 
  unnest(GO_term)
new <- new[!new$GO_term=="-",]

write.table(new, file="louse_GO_annotations.txt", sep="\t", quote = F, row.names = F,
            col.names = T)
