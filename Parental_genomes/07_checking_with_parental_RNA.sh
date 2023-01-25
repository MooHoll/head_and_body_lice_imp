#---------------------------------------------------------------------------------------
# Check to see if the called unique SNPs show >0.9 read propoertion in the pure parental RNA-seq
#---------------------------------------------------------------------------------------

module load bedtools/2.28.0
module load samtools/1.9

# Count original number of aligned reads
samtools view -c -F 260 BB_M_1.STAR.genome.bam  # 106643158
samtools view -c -F 260 BB_M_2.STAR.genome.bam  # 103815440
samtools view -c -F 260 HH_M_1.STAR.genome.bam  # 108379640
samtools view -c -F 260 HH_M_2.STAR.genome.bam  # 143988726

# Need to unzip all vcfs for this to work, urgh
module load htslib/1.9
bgzip -d snp_for_masking.vcf.gz 
grep -v "#" snp_for_masking.vcf | wc -l # 360525

# Make a snp region file
module load bedops/2.4.26 
vcf2bed < snp_for_masking.vcf > all_snps.bed
paste <(cut -f1 all_snps.bed) <(cut -f2 all_snps.bed) <(cut -f2 all_snps.bed) > snp_regions.bed

# Sort and index bams
samtools sort BB_M_1.STAR.genome.bam > sorted_BB_1.bam
samtools sort BB_M_2.STAR.genome.bam > sorted_BB_2.bam
samtools sort HH_M_1.STAR.genome.bam > sorted_HH_1.bam
samtools sort HH_M_2.STAR.genome.bam > sorted_HH_2.bam

samtools index sorted_BB_1.bam
samtools index sorted_BB_2.bam
samtools index sorted_HH_1.bam
samtools index sorted_HH_2.bam

# Pull out actual called SNPs for later checking
cut -f1,2,7 snp_regions.bed > all_snps.txt

#----------------------
# Count nucleotides at each position
#----------------------
#!/bin/bash

#PBS -N counting_snps_in_RNAseq
#PBS -l walltime=20:00:00
#PBS -l vmem=50gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=1

# Run script in the working directory it was submitted in  
cd $PBS_O_WORKDIR 

# Load modules
module load samtools/1.9
module load bam-readcount/20161006 

for file in $(ls *bam)
do
	base=$(basename ${file} ".bam")
    bam-readcount -l snp_regions.bed -f JCVI_LOUSE_1.0_genomic.fa ${file} > ${base}_snp_counts.txt
done