# Make an N masked genome for each colony

module load bedtools/2.28.0
module load bcftools/1.9
module load htslib/1.9 
module load samtools/1.9

# Need one SNP file which masks the unique SNPs of both ecotypes to avoid mapping bias
bgzip bl_unique.vcf
bgzip hl_unique.vcf

# Need index for the merge command
bcftools index bl_unique.vcf.gz
bcftools index hl_unique.vcf.gz

# Make one SNP file per colony
bcftools merge -o snp_for_masking.vcf.gz bl_unique.vcf.gz hl_unique.vcf.gz

# Make an N masked genome per colony
bedtools maskfasta -fi JCVI_LOUSE_1.0_genomic.fa -bed snp_for_masking.vcf.gz -fo louse_masked_genome.fa 

# Note: can use "snp_for_masking.vcf.gz" for pulling out info with bam read count, using all available positions then