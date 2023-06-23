#!/bin/bash

#PBS -N alignment
#PBS -l walltime=03:00:00
#PBS -l vmem=10gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=16

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load software needed 
module load bowtie2/2.3.5.1
module load samtools/1.8

echo "making the genome index"
bowtie2-build /scratch/monoallelic/hjm32/lice/genome/JCVI_LOUSE_1.0_genomic.fa louse

echo "starting alignment"
for file in $(ls *_1.fq.gz)
do
	base=$(basename $file ".f0.trimmed_1.fq.gz")
    bowtie2 --sensitive --threads 14 \
    -x louse \
    -1 ${base}.f0.trimmed_1.fq.gz -2 ${base}.f0.trimmed_2.fq.gz \
    -S ${base}.sam
done

echo "convert sams to bams"
for file in $(ls *.sam)
do
	base=$(basename $file ".sam")
    samtools view -bS ${base}.sam > ${base}.bam
done

echo "sort bams"
for file in $(ls *.bam)
do
	base=$(basename $file ".bam")
    samtools sort -@ 14 -o ${base}_sorted.bam ${file}
done