#!/bin/bash -l

#SBATCH --job-name=minimap2
#SBATCH -D /group/jrigrp6/andropogon_sequence
#SBATCH -e /group/jrigrp6/andropogon_sequence/slurm_log/minimap2_sterror_%j.txt
#SBATCH -o /group/jrigrp6/andropogon_sequence/slurm_log/minimap2_stdoutput_%j.txt
#SBATCH -p med2
#SBATCH -t 7-00:00
#SBATCH --cpus-per-task 8


module load minimap2

# Sorghum vs Andro
#minimap2 -cx asm10 --cs -t 8 \
#~/Andropogoneae/sorghum_v_andropogon/GCA_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna \
#data/genome/corteva_andro/Andropogon_geardii_hifiasm-bionano-nohets_v2.fasta \
#> data/minimap/andro_sorg_alignment.paf


# Andro vs Andro
minimap2 -cx asm5 --cs -t 8 \
data/genome/corteva_andro/Andropogon_geardii_hifiasm-bionano-nohets_v2.fasta \
data/genome/corteva_andro/Andropogon_geardii_hifiasm-bionano-nohets_v2.fasta \
> data/minimap/andro_andro_alignment.paf
