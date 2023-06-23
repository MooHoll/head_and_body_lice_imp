#!/bin/bash

#PBS -N trimming
#PBS -l walltime=05:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=10

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load module
module load cutadapt/1.11
module load pigz/2.3.3

for file in $(ls *_1.fq.gz)
do
	base=$(basename $file "_1.fq.gz")
	cutadapt \
	-u 10 -U 10 \
	-o ${base}_trimmed_1.fq.gz \
	-p ${base}_trimmed_2.fq.gz \
	${base}_1.fq.gz \
	${base}_2.fq.gz
done