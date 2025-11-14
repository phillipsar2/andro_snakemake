#!/bin/bash -l

#SBATCH --job-name=fastqc
#SBATCH -D /group/jrigrp6/andropogon_sequence
#SBATCH -e /group/jrigrp6/andropogon_sequence/slurm_log/sterror_fastqc_%j.txt
#SBATCH -o /group/jrigrp6/andropogon_sequence/slurm_log/stdoutput_fastqc_%j.txt
#SBATCH -p high2
#SBATCH -t 30-00:00
#SBATCH --mem 32G


module load bio3

#fastqc -o qc/fastqc -f fastq /group/jrigrp10/andropogon_shortreads/ucd_seq/*.fastq.gz
#fastqc -o qc/fastqc -f fastq /group/jrigrp10/andropogon_shortreads/IN*.fastq.gz
fastqc -o qc/fastqc -f fastq data/raw/trimmed/*.fastq.gz
