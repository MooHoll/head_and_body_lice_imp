# ---------------------------------------------
# Creating psudoreplicates 
# ---------------------------------------------

# Based off this paper:
# https://doi.org/10.1038/s41598-022-11302-9

# Following this link to select reads from fastq files with replacement
# creating a new file of equal size, as a psuedoreplicate
# https://thelegendofbioinformatics.wordpress.com/2017/08/24/random-fastq-subsampling/

module load python/gcc/3.5.5

#  -s = make new files 99% of the size of the original (100% doesn't want to work)
#  -n = male this many new paired ended replicates
#  -r = with replacement (this part is the most important to generate some variation)
python randomReadSubSample.py -f1 SRR9617639_trim_1.fq -f2 SRR9617639_trim_2.fq -s 0.99 -n 2 -r 1 -o SRR9617639_replicate
python randomReadSubSample.py -f1 SRR9617640_trim_1.fq -f2 SRR9617640_trim_2.fq -s 0.99 -n 2 -r 1 -o SRR9617640_replicate