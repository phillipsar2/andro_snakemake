#!/bin/bash -l

#SBATCH --job-name=minimap2
#SBATCH -D /group/jrigrp6/andropogon_sequence
#SBATCH -e /group/jrigrp6/andropogon_sequence/slurm_log/minimap2_sterror_%j.txt
#SBATCH -o /group/jrigrp6/andropogon_sequence/slurm_log/minimap2_stdoutput_%j.txt
#SBATCH -p bigmemm
#SBATCH -t 7-00:00
#SBATCH --cpus-per-task 8


module load minimap2

# Sorghum vs Andro
minimap2 -cx asm10 --cs -t 8 \
~/Andropogoneae/sorghum_v_andropogon/GCA_000003195.3_Sorghum_bicolor_NCBIv3_genomic.fna \
data/genome/JGIgenome/v1/Andropogon_gerardii_var_Kellogg_1272_HAP1_V1_release/Andropogon_gerardii_var_Kellogg_1272/sequences/Andropogon_gerardii_var_Kellogg_1272.mainGenome.fasta \
> data/minimap/JGIandro_sorg_alignment_20211015.paf


# Andro vs Andro
#minimap2 -cx asm5 --cs -t 8 \
#data/genome/corteva_andro/Andropogon_geardii_hifiasm-bionano-nohets_v2.fasta \
#data/genome/corteva_andro/Andropogon_geardii_hifiasm-bionano-nohets_v2.fasta \
#> data/minimap/andro_andro_alignment.paf
