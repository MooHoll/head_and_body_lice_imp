#!/bin/bash

#PBS -N trimming
#PBS -l walltime=01:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=10

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load module
module load cutadapt/1.11
module load pigz/2.3.3

for file in $(ls *_1.fastq)
do
	base=$(basename $file "_1.fastq")
	cutadapt \
	-u 10 -U 10 \
	-o ${base}_trim_1.fq \
	-p ${base}_trim_2.fq \
	${base}_1.fastq \
	${base}_2.fastq
done