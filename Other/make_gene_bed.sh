# Make gene bed file from genome annotation

awk '$3 == "gene"' JCVI_LOUSE_1.0_genomic.gff > genes
cut -f1,4,5,9 genes > louse_genes.bed
sed -i 's/;.*//g' louse_genes.bed
sed -i 's/ID=//g' louse_genes.bed