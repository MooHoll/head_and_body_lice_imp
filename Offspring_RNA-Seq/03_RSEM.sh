#!/bin/bash

#PBS -N rsem
#PBS -l walltime=12:00:00
#PBS -l vmem=50gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=16

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load module
module load rsem/1.3.0
module load star/2.7.9a

# Renamed files to remove trimmed word
# for file in *.gz ; do mv "${file}" "${file/_trimmed/}"; done

#echo "convert gff to gtf for RSEM"
#module load cufflinks/2.2.1
#gffread JCVI_LOUSE_1.0_genomic.gff -T -o JCVI_LOUSE_1.0_genomic.gtf

#echo "prepapre genome"
#rsem-prepare-reference --gtf JCVI_LOUSE_1.0_genomic.gtf --star -p 1 louse_masked_genome.fa louse

echo "run RSEM"
for file in $(ls *_1.fq.gz)
do
    base=$(basename ${file} "_1.fq.gz")
    rsem-calculate-expression --star --star-gzipped-read-file -p 16 --star-output-genome-bam \
    --paired-end ${base}_1.fq.gz ${base}_2.fq.gz \
    /scratch/monoallelic/hjm32/lice/genome/louse ${base}
done
