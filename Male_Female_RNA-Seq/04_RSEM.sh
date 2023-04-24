#!/bin/bash

#PBS -N rsem
#PBS -l walltime=02:00:00
#PBS -l vmem=50gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=16

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load module
module load rsem/1.3.0
module load star/2.7.9a

#echo "convert gff to gtf for RSEM"
#module load cufflinks/2.2.1
#gffread JCVI_LOUSE_1.0_genomic.gff -T -o JCVI_LOUSE_1.0_genomic.gtf

#echo "prepapre genome"
#rsem-prepare-reference --gtf JCVI_LOUSE_1.0_genomic.gtf --star -p 1 JCVI_LOUSE_1.0_genomic.fa louse

echo "run RSEM"
for file in $(ls *1.fq)
do
    base=$(basename ${file} "1.fq")
    rsem-calculate-expression --star -p 16 --star-output-genome-bam \
    --paired-end ${base}1.fq ${base}2.fq \
    /scratch/monoallelic/hjm32/lice/genome/louse ${base}
done
