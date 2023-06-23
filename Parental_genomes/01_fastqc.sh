#!/bin/bash

#PBS -N fastqc
#PBS -l walltime=01:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

cd $PBS_O_WORKDIR

module load fastqc/0.11.5

# Create a list of the files to be called
for file in $(ls *.gz)
do
  	fastqc -t 7 ${file}
done