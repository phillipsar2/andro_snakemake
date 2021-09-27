#!/bin/bash -l

#SBATCH --job-name=jellyfish
#SBATCH -D /group/jrigrp6/andropogon_sequence
#SBATCH -e /group/jrigrp6/andropogon_sequence/slurm_log/jellyfish_sterror_%j.txt
#SBATCH -o /group/jrigrp6/andropogon_sequence/slurm_log/jellyfish_stdoutput_%j.txt
#SBATCH -p bigmemh
#SBATCH -t 7-00:00
#SBATCH --cpus-per-task 10


# Alyssa Phillips
# 8/3/21
# Compute k-mer frequencies with Jellyfish


module load jellyfish/2.2.6

#sample = /group/jrigrp10/andropogon_shortreads/INDY_PCRfree_Loretta_BB-2

sample=/group/jrigrp10/andropogon_shortreads/INDQ_PCRfree_Loretta_FUL-3_33_AAGATACT_Andropogon_gerardii.merge

r1="$sample".R1.fastq.gz
r2="$sample".R2.fastq.gz


# Count k-mers
# -C : canonical k-mers
# -m : k-mer length - 21 recommended by GenomeScope
# -s : memory (currently set 1GB)
# -t : threads
# -F : number of files to read simultaneously
#jellyfish count -C -m 21 -F 2 -s 1000000000 -t 10 <(zcat /group/jrigrp10/andropogon_shortreads/INDY_PCRfree_Loretta_BB-2_19_CATAATAC_Andropogon_gerardii.merge.R1.fastq.gz) <(zcat /group/jrigrp10/andropogon_shortreads/INDY_PCRfree_Loretta_BB-2_19_CATAATAC_Andropogon_gerardii.merge.R2.fastq.gz) -o data/jellyfish/INDY_reads.jf
jellyfish count -C -m 21 -F 2 -s 1000000000 -t 10 $r1 $r2 -o data/jellyfish/INDQ_reads.jf


# Export k-mer count histogram
jellyfish histo -t 10 data/jellyfish/INDQ_reads.jf > data/jellyfish/INDQ_reads.histo
