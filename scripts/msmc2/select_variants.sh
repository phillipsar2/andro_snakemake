#!/bin/bash
#SBATCH --partition=high2
#SBATCH --job-name=make_join_gvcf
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb
#SBATCH --time=150:00:00
#SBATCH --output=gvcf_%A_%a.log
#SBATCH --array 7

module load GATK/4.2.3.0
module load samtools

cd /group/jrigrp11/nam_gvcfs

gatk SelectVariants \
    -R Zm-B73-REFERENCE-NAM-5.0.fa \
    -V gendb://nam_db_"$SLURM_ARRAY_TASK_ID" \
    -L $SLURM_ARRAY_TASK_ID:1-70000000 \
    -L $SLURM_ARRAY_TASK_ID:90000000-180000000 \
    -O nam_"$SLURM_ARRAY_TASK_ID".vcf
