#!/bin/bash -l

#SBATCH --job-name=pcangsd
#SBATCH -D /group/jrigrp6/andropogon_sequence
#SBATCH -e /group/jrigrp6/andropogon_sequence/slurm_log/sterror_pca_%j.txt
#SBATCH -o /group/jrigrp6/andropogon_sequence/slurm_log/stdoutput_pca_%j.txt
#SBATCH -p high2
#SBATCH -t 7-00:00
#SBATCH --mem 256G

module load angsd
module load bio3
module load pcangsd

source activate pcangsd

pcangsd -beagle data/angsd/pca/highcov.merged.8dp70.beagle.gz -threads 30 -minMaf 0.2 -o data/angsd/pca/highcov.merged.8dp70
