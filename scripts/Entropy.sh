#!/bin/bash -l

#SBATCH --job-name=entropy
#SBATCH -D /group/jrigrp6/andropogon_sequence
#SBATCH -e /group/jrigrp6/andropogon_sequence/slurm_log/sterror_entropy%a_%j.txt
#SBATCH -o /group/jrigrp6/andropogon_sequence/slurm_log/stdoutput_entropy%a_%j.txt
#SBATCH -p med2
#SBATCH -t 30-00:00
#SBATCH --mem 8G
#SBATCH --array=2-14

module load bio3
source activate entropy-2.0

pre=data/entropy/lowcov/10k.lowcov.miss5.edit
file=10k.lowcov.miss5.edit

ngroups=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 )

chain=c1
k=${ngroups[$SLURM_ARRAY_TASK_ID]}

out=/scratch/aphillip/entropy/"$file".k"$k"."$chain".hdf5

echo $out

module list

mkdir -p  /scratch/aphillip/entropy

# -D 1 calculate WAIC estimates
# -t store every 5th step after burn-in
# -l 50000 run 50k MCMC steps
# -r set seed

# create empty file
touch $out

# run
entropy -i "$pre".mpgl \
-n "$pre".ploidy_inds.txt \
-k $k \
-q "$pre".qk"$k"inds.txt \
-t 5 \
-D 1 \
-l 70000 \
-r 22 \
-o $out

mv $out data/entropy/lowcov/"$file".k"$k"."$chain".hdf5
#rm -r mkdir /scratch/aphillip/entropy

# grab WAIC
estpost.entropy -s 3 -p likdata data/entropy/lowcov/"$file".k"$k"."$chain".hdf5

# grab q
estpost.entropy -p q -s 0 data/entropy/lowcov/"$file".k"$k"."$chain".hdf5 -o data/entropy/lowcov/"$file"."$pre".qest."$k"."$chain".txt

