#---------------------------------------------------------------------------------------
# Count reads per informative SNP in the offspring from crosses
#---------------------------------------------------------------------------------------

#!/bin/bash

#PBS -N counting_reads
#PBS -l walltime=30:00:00
#PBS -l vmem=50gb
#PBS -m bea
#PBS -M hollie_marshall@hotmail.co.uk
#PBS -l nodes=1:ppn=1

# Run script in the working directory it was submitted in  
cd $PBS_O_WORKDIR 

# Load modules
module load samtools/1.9
module load bam-readcount/20161006 

# Format the bams
for file in $(ls *bam)
do
    base=$(basename ${file} ".STAR.genome.bam")
    echo "counting reads"
    samtools view -c -F 260 ${file}
    echo "sorting"
    samtools sort ${file} > ${base}_sorted.bam
    echo "indexing"
    samtools index ${base}_sorted.bam
done

# bh_1 113066526
# bh_4 108159116
# bh_5 129713200
# bh_6 121473822
# hb_1 124896206
# hb_2 129068584
# hb_3 128443338
# hb_5 135815786
# hb_6 108814356

# Count the snps
echo "counting SNPs"
for file in $(ls *_sorted.bam)
do
	base=$(basename ${file} "_sorted.bam")
    bam-readcount -l snp_regions.bed -f JCVI_LOUSE_1.0_genomic.fa ${file} > ${base}_snp_counts.txt
done