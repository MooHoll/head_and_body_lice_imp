#!/bin/bash

#PBS -N snp_calling
#PBS -l walltime=00:05:00
#PBS -l vmem=50gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=1
#PBS -q devel

# Run script in the working directory it was submitted in  
cd $PBS_O_WORKDIR 

# Load software needed
module load samtools/1.9 
module load vcftools/0.1.14
module load freebayes/1.1.0

# NOTE: freebayes is single thread, can multithread with freebayes-parallel
# but this only splits the genome using a simple python script, will just avoid

echo "indexing bams"
for file in $(ls *_sorted_RG_deduplicated.bam)
do
    samtools index ${file}
done

echo "call SNPs"
# min count 2 reads for alternative alleles, 
# min 5 reads per SNP, ignore complex events, indels and mnps
for file in $(ls *_sorted_RG_deduplicated.bam)
do
  	base=$(basename ${file} "_sorted_RG_deduplicated.bam")
    freebayes \
        -f JCVI_LOUSE_1.0_genomic.fa \
        -C 2 \
        --min-coverage 5 \
        -u \
        -i \
        -X \
        -b ${file} \
        > ${base}_vcf.gz
done

echo "filtering SNPs"
for file in $(ls *gz)
do
	base=$(basename ${file} "vcf.gz")
    vcftools --gzvcf ${file} --max-alleles 2 --minQ 20 --min-meanDP 10 \
    --recode --recode-INFO-all --out ${base}
done