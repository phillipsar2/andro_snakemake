#!/bin/bash -l

#SBATCH --job-name=smudgeplot
#SBATCH -D /group/jrigrp6/andropogon_sequence
#SBATCH -e /group/jrigrp6/andropogon_sequence/slurm_log/smudgeplot_sterror_%j.txt
#SBATCH -o /group/jrigrp6/andropogon_sequence/slurm_log/smudgeplot_stdoutput_%j.txt
#SBATCH -p bigmemh
#SBATCH -t 7-00:00
#SBATCH --mem 200G
#SBATCH --cpus-per-task 12

#  https://github.com/KamilSJaron/smudgeplot

sample=INDE_PCRfree_Loretta_SAL-5_6_GCACGGAC_Andropogon_gerardii

## Calculate kmer frequences from the fasta files and then store in a histogram
## Takes ~40 min to get the histogram

module load KMC

#mkdir /scratch/aphillip/"$sample"

#ls /group/jrigrp10/andropogon_shortreads/"$sample"*.fastq.gz > data/kmer/"$sample".FILES

#kmc -k21 -t12 -m64 -ci1 -cs10000 \
#@data/kmer/"$sample".FILES \
#data/kmer/"$sample" \
#/scratch/aphillip/"$sample"

#kmc_tools transform data/kmer/"$sample" \
#histogram data/kmer/"$sample".hist -cx10000

#rm -r /scratch/aphillip/"$sample"

## extract genomic kmers within a coverage threshold then generate smudgeplot

module load smudgeplot
source activate smudgeplot-0.2.4

#L=$(smudgeplot.py cutoff data/kmer/"$sample".hist L)
#U=$(smudgeplot.py cutoff data/kmer/"$sample".hist U)
#echo $L $U > data/kmer/"$sample".limits # these need to be sane values
        # L should be like 20 - 200
        # U should be like 500 - 3000

## first tried 20 and 1900, really large file (23.2 G)
## 50, 1500 is a much smaller file (6.6G)

L=50
U=1500

kmc_tools transform data/kmer/"$sample" -ci"$L" -cx"$U" \
dump -s data/kmer/"$sample"_L"$L"_U"$U".dump


## hetkmers requires memory 15x the size of the .dump file -- needed 200G, 3hr run time
## can increase L and lower U in order to decrease the file size
smudgeplot.py hetkmers -o data/kmer/"$sample"_L"$L"_U"$U" < \
data/kmer/"$sample"_L"$L"_U"$U".dump

## memory requirement = 1G at most
smudgeplot.py plot data/kmer/"$sample"_L"$L"_U"$U"_coverages.tsv \
-o data/smudgeplot/"$sample"

source deactivate

