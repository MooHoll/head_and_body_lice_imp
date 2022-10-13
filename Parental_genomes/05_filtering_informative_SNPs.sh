# -----------------------------------------------
# Get homozygous alternative SNPS
# -----------------------------------------------

# Count total heterozygote SNPs 
gunzip *recode.vcf.gz
grep -c 0/1 <sample.vcf>
# bl 484226
# hl 358891

# Count total homozygous SNPs
grep -c 1/1 <sample.vcf>
# bl 306776
# 355073

# Write file containing only the alternate homozygous SNPs 
grep -e ^# -e "1/1" bl.recode.vcf > bl_alt_homozygous_snps.vcf
grep -e ^# -e "1/1" hl.recode.vcf > hl_alt_homozygous_snps.vcf

# Make files with SNPs unique to each ecotype
module load bedtools/2.28.0

# Remove uninformative SNPs that males and queens have in common
# -A needed to remove whole SNP information if there is overlap, -header needed for future commands to recognise the .vcf format 
subtractBed -header -A -a bl_alt_homozygous_snps.vcf -b hl_alt_homozygous_snps.vcf > bl_unique.vcf
grep -v ^# bl_unique.vcf | wc -l
# 156114

subtractBed -header -A -a hl_alt_homozygous_snps.vcf -b bl_alt_homozygous_snps.vcf > hl_unique.vcf
grep -v ^# hl_unique.vcf | wc -l
# 204411

