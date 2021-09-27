#!/bin/sh

#SBATCH --job-name=kmc
#SBATCH -D /group/jrigrp6/andropogon_sequence
#SBATCH -e /group/jrigrp6/andropogon_sequence/slurm_log/kmc_sterror_%j.txt
#SBATCH -o /group/jrigrp6/andropogon_sequence/slurm_log/kmc_stdoutput_%j.txt
#SBATCH -p bigmemh
#SBATCH -t 7-00:00
#SBATCH --cpus-per-task 16
#SBATCH --mem 64G


# Alyssa Phillips
# 8/11/21
# Compute k-mer frequencies with KMC for use in Smudgeplot

module load KMC/3.1.1

sample=/group/jrigrp10/andropogon_shortreads/INDQ_PCRfree_Loretta_FUL-3_33_AAGATACT_Andropogon_gerardii.merge

r1="$sample".R1.fastq.gz
r2="$sample".R2.fastq.gz

# database name
db=data/kmc/INDQ_PCRfree_Loretta_FUL-3_33_AAGATACT_Andropogon_gerardii.merge
tmp=/scratch/aphillip/kmc
hist=data/kmc/INDQ_PCRfree_Loretta_FUL-3_33_AAGATACT_Andropogon_gerardii.merge_k21.hist

#mkdir /scratch/aphillip/kmc

#ls "$sample"* > FILES

# (1) calculate kmer frequencies
# kmer 21, 16 threads, 64G of memory, counting kmer coverages between 1 and 10000x
# -ci : excludes kmers occuring less than this many times
# -cs : maximum counter value
#kmc -k21 -t16 -m64 -ci1 -cs10000 @FILES $db $tmp

# (2)  generate histogram of kmer frequencies
# -cx : maximum value of counter to be stored in the histogram file
#kmc_tools transform $db histogram $hist -cx10000

# (3) evaluate upper and lower kmer bounds - manually with smudgeplot
L=15
U=3300


# (4) extract kmers in coverage range
#kmc_tools transform $db -ci"$L" -cx"$U" reduce "$db"_L"$L"_U"$U"
kmc_tools transform $db -ci"$L" -cx"$U" dump -s "$db"_L"$L"_U"$U".dump


# (5) compute set of kmer pairs for reduced file
