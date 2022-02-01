#!/bin/bash -l

#SBATCH --job-name=ebg
#SBATCH -D /group/jrigrp6/andropogon_sequence
#SBATCH -e /group/jrigrp6/andropogon_sequence/slurm_log/sterror_ebg%j.txt
#SBATCH -o /group/jrigrp6/andropogon_sequence/slurm_log/stdoutput_ebg%j.txt
#SBATCH -p bigmemh
#SBATCH -t 7-00:00
#SBATCH --mem 2G

## Run the diseq model in EBG to caluclate genotype likelihoods ##

sample=Chr07B.9x.10k
sample_count=14
ploidy=3
loci=9290

total=data/ebg/total_reads."$sample".txt
alt=data/ebg/alt_reads."$sample".txt
error=data/ebg/error_rate."$sample".txt

prefix=data/ebg/"$sample"

~/toolsfordayz/polyploid-genotyping/ebg/ebg diseq \
-t $total \
-a $alt \
-e $error \
-p $ploidy \
-n $sample_count \
-l $loci \
--iters 1000 \
--prefix $prefix

# paper states it typically converges in 50-100 interations
