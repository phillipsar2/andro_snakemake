import config
import pandas as pd

# Different SAMPLE names
SAMPLE = glob_wildcards("/group/jrigrp10/andropogon_shortreads/{sample}.merge.R1.fastq.gz").sample
#SAMPLE = glob_wildcards("/group/jrigrp10/andropogon_shortreads/{sample}_1.fq.gz").sample
#print(SAMPLE)

# BAM names encompass all samples regardless of their fastq pattern
BAM = glob_wildcards("data/interm/mark_dups/{bam}.dedup.bam").bam
#print(BAM)

# MERGE contains a list of the bams that belong to each genotype (GENO)  so they can be merged
file = pd.read_csv("bams_to_merge.tsv", sep = "\t", header = 0)
MERGE_A = list(file.Merge_A)
MERGE_B = list(file.Merge_B)
GENO = list(file.Genotype)

#print(MERGE_A)
#print(MERGE_B)
#print(GENO)

# TEST contains the samples for which AB should be calculate to determine ploidy
# 9x, 6x, unknown
#TEST = ['INDE_PCRfree_Loretta_SAL-5_6_GCACGGAC_Andropogon_gerardii', 'INDN_PCRfree_5-3-McKain_254-WI-Diversity_46_ATATGGAT_Andropogon_gerardii', 'INDA_PCRfree_Texas-2_41_CCAAGTCT_Andropogon_gerardii']
filep = pd.read_csv("bams_nquire.csv", sep = "\t", header = None)
UNKNOWN = list(filep[0])


# List of genotypes and ploidy levels
# BAM = 
# ploidy = 

# Number of intervals for GATK SplitIntervals
INTERVALS = ["{:04d}".format(x) for x in list(range(200))]


# Rule all describes the final output of the pipeline
rule all:
    input:
        ## Aligning reads
#        expand("data/interm/mark_dups/{sample}.dedup.bam", sample = SAMPLE),
#        expand("reports/bamqc/{bam}_stats/qualimapReport.html", bam = BAM),
        ## Merging low-coverage bams
#        expand("data/interm/mark_dups/{geno}.{merge_A}.{merge_B}.merged.dedup.bam", zip, merge_A = MERGE_A, merge_B = MERGE_B, geno = GENO),
#        expand("data/interm/mark_dups/{geno}.{merge_A}.{merge_B}.merged.rg.dedup.bam", zip, merge_A = MERGE_A, merge_B = MERGE_B, geno = GENO),
        ## nQuire
#        hist = expand("data/nQuire/{unknown}_denoised.hist", unknown = UNKNOWN),
        r = "data/nQuire/nquire_results_denoised.txt",
        ## Calculate AB
#        expand("data/ab/vcf/{test}.vcf", test = TEST),
#        expand(""data/ab/table/{test}.AD.table", test = TEST),
        ## SNP Calling
#        expand("data/gvcf/{bam}.p{ploidy}.g.vcf", zip, ploidy = PLOIDY, bam = BAM)

# Rules
include: "rules/mapping.smk"
#include: "rules/process_bam.smk"
include: "rules/determine_ploidy.smk"
#include: "rules/calc_AB.smk"
#include: "rules/calling_snps.smk"
