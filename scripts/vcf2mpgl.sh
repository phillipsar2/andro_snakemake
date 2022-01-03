#!/bin/bash -l

#SBATCH --job-name=vcf2mpgl
#SBATCH -D /group/jrigrp6/andropogon_sequence
#SBATCH -e /group/jrigrp6/andropogon_sequence/slurm_log/sterror_vcf2mpgl_%j.txt
#SBATCH -o /group/jrigrp6/andropogon_sequence/slurm_log/stdoutput_vcf2mpgl_%j.txt
#SBATCH -p high2
#SBATCH -t 30-00:00
#SBATCH --mem 4G

module load bio3 # load R packages
module load spack/R/4.1.0 # load vcfR

source activate entropy-2.0 # enter entropy environment

# input vcf file
vcf=final.1000ln.test.vcf.gz
#vcf=/group/jrigrp10/tripsacum_dact/data/vcf/final/merged/final.1000ln.test.vcf.gz

# Covert vcf to mpgl
Rscript scripts/vcf2mpgl.R $vcf
