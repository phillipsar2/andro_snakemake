#!/bin/bash -l

#SBATCH --job-name=bwa_mem
#SBATCH -D /group/jrigrp6/andropogon_sequence
#SBATCH -e /group/jrigrp6/andropogon_sequence/slurm_log/sterror_mapping_%j.txt
#SBATCH -o /group/jrigrp6/andropogon_sequence/slurm_log/stdoutput_mapping_%j.txt
#SBATCH -p bigmemh
#SBATCH -t 30-00:00
#SBATCH --mem 32G
#SBATCH --cpus-per-task 8

module load bwa-mem2/2.2.1_x64-linux

bwa-mem2 mem -t 8 data/genome/corteva_andro/Andropogon_geardii_hifiasm-bionano-nohets_v2.fasta \
/group/jrigrp10/andropogon_shortreads/AN20N108*_1.fq.gz \
/group/jrigrp10/andropogon_shortreads/AN20N108*_2.fq.gz | \
samtools view -Sb > data/interm/mapped_bam/AN20N108_CKDL200149283-1a-7UDI297-AK22043_HNLFVDSXX_L1.mapped.bam

