#!/bin/bash -l

#SBATCH --job-name=bwa_index
#SBATCH -D /group/jrigrp6/andropogon_sequence
#SBATCH -e /group/jrigrp6/andropogon_sequence/slurm_log/sterror_%j.txt
#SBATCH -o /group/jrigrp6/andropogon_sequence/slurm_log/stdoutput_%j.txt
#SBATCH -p bigmemm
#SBATCH -t 7-00:00
#SBATCH --mem 124G

module load bwa-mem2/2.2.1_x64-linux

bwa-mem2 index data/genome/JGIgenome/v1/Andropogon_gerardii_var_Kellogg_1272_HAP1_V1_release/Andropogon_gerardii_var_Kellogg_1272/sequences/Andropogon_gerardii_var_Kellogg_1272.mainGenome.fasta
