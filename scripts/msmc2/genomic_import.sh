#SBATCH --partition=high2
#SBATCH --job-name=dbimport
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=40gb
#SBATCH --time=150:00:00
#SBATCH --output=serial_test_%A_%a.log
#SBATCH --array 1-10

module load GATK/4.2.3.0

cd /group/jrigrp10/ager

gatk --java-options "-Xmx40g -Xms8g" GenomicsDBImport \
-V 1B.gvcf \
-V 2B.gvcf \
-V 1C.gvcf \
-V 2C.gvcf \
-V 2A.gvcf \
--genomicsdb-workspace-path ager_db_"$SLURM_ARRAY_TASK_ID"  \
--genomicsdb-vcf-buffer-size 131072000 \
--overwrite-existing-genomicsdb-workspace TRUE \
-L $SLURM_ARRAY_TASK_ID
