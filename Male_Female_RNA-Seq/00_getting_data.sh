#!/bin/bash

#PBS -N getting_data
#PBS -l walltime=01:00:00
#PBS -l vmem=20gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=8

# Change directory to the one the job was submitted in
cd $PBS_O_WORKDIR 

# Load required modules
module load sratoolkit/2.11.1

# Download data, from paper: Investigation into sex differentiation pathway of hemimetabolous insects
# female
prefetch SRR9617639
# male
prefetch SRR9617640

# Move all files out of the directories

# Make SRR files into fastq
fasterq-dump --split-files --split-3 SRR9617639.sra
fasterq-dump --split-files --split-3 SRR9617640.sra