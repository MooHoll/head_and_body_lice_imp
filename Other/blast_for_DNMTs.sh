#-------------------------------------------------------------------
# Reciprocal blast to identify the DNMT1 and 3 genes in the louse
#-------------------------------------------------------------------

# Get protein seqs of all louse genes
conda create -n agat_env
conda activate agat_env
conda install -c bioconda agat

agat_sp_extract_sequences.pl --gff JCVI_LOUSE_1.0_genomic.gff -f JCVI_LOUSE_1.0_genomic.fa -p -o louse_proteins.fa

# Make a blast databases
module load blast+/2.9.0

makeblastdb -in louse_proteins.fa -parse_seqids -dbtype prot

# Using http://v2.insect-genome.com/Pcg, could download the protein .fa seq in one file
# for DNMT1 and DNMT3a for 321 and 110 species respectively 
makeblastdb -in dnmt1_proteins_321_insect_species.fa -parse_seqids -dbtype prot
makeblastdb -in dnmt3a_proteins_110_insect_species.fa -parse_seqids -dbtype prot

# reciprocal for DNMT1
blastp -query louse_proteins.fa \
-db "dnmt1_proteins_321_insect_species.fa" -evalue 1e-3 -max_target_seqs 1 \
-outfmt 6 -out louse_DNMT1.txt

blastp -query dnmt1_proteins_321_insect_species.fa \
-db "louse_proteins.fa" -evalue 1e-3 -max_target_seqs 1 \
-outfmt 6 -out louse_DNMT1_reciprocol.txt

module load R/4.0.0
R
library(readr)
data <- read_delim("louse_DNMT1_reciprocol.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names=F)
table(data$X2)
#rna-XM_002431833.1 (gene-Phum_PHUM556380) rna-XM_002432115.1 (gene-Phum_PHUM574160)
#               343                 17 


# reciprocal for DNMT3
blastp -query louse_proteins.fa \
-db "dnmt3a_proteins_110_insect_species.fa" -evalue 1e-3 -max_target_seqs 1 \
-outfmt 6 -out louse_DNMT3.txt

blastp -query dnmt3a_proteins_110_insect_species.fa \
-db "louse_proteins.fa" -evalue 1e-3 -max_target_seqs 1 \
-outfmt 6 -out louse_DNMT3_reciprocol.txt

R
library(readr)
data <- read_delim("louse_DNMT3_reciprocol.txt", "\t", escape_double = FALSE, trim_ws = TRUE, col_names=F)
table(data$X2)
# rna-XM_002423125.1 (gene-Phum_PHUM041250) rna-XM_002425349.1 rna-XM_002425839.1 rna-XM_002426513.1 
#                50                  1                  1                  3 
#rna-XM_002426548.1 rna-XM_002428121.1 rna-XM_002429418.1 rna-XM_002429731.1 
#                 1                  4                  2                  3 
#rna-XM_002429883.1 (gene-Phum_PHUM456660) rna-XM_002431833.1 
#                15                  3 