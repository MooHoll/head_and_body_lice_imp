#!/bin/bash

#PBS -N fastqc
#PBS -l walltime=01:00:00
#PBS -l vmem=10gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=10

# Run script in the working directory it was submitted in 
cd $PBS_O_WORKDIR 

# Load module
module load fastqc/0.11.5

for file in $(ls *.fastq)
do
	fastqc ${file}
done
