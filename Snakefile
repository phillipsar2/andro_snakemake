import config
import pandas as pd

# Different SAMPLE names
SAMPLE = glob_wildcards("/group/jrigrp10/andropogon_shortreads/{sample}.merge.R1.fastq.gz").sample
#SAMPLE = glob_wildcards("/group/jrigrp10/andropogon_shortreads/{sample}_1.fq.gz").sample
#print(SAMPLE)

# BAM names encompass all samples regardless of their fastq pattern
#BAM = glob_wildcards("data/interm/mark_dups/{bam}.dedup.bam").bam
#BAM = glob_wildcards("data/final_bams/lowcov/{bam}.dedup.bam").bam
BAM = glob_wildcards("data/final_bams/highcov/{bam}.dedup.bam").bam
#print(BAM)

# MERGE contains a list of the bams that belong to each genotype (GENO)  so they can be merged
file = pd.read_csv("bams_to_merge.tsv", sep = "\t", header = 0)
MERGE_A = list(file.Merge_A)
MERGE_B = list(file.Merge_B)
GENO = list(file.Genotype)

# HIGH contains the bams to be subsampled to low coverage
HIGH = glob_wildcards("data/interm/mark_dups/IN{high}.dedup.bam").high

# filep contains the bams nQuire was run on to determine ploidy
filep = pd.read_csv("bams_nquire.csv", sep = "\t", header = None)
UNKNOWN = list(filep[0])

# List of genotypes and ploidy levels
fileg = pd.read_csv(config.highcov, header = 0)
PLOIDY = list(fileg.Ploidy)
VCF = list(fileg.Sequencefile)

# List of chromsomes
chr = pd.read_csv(config.contig_list, header = None)
CHROM = list(chr[0])

# Select working with high or low coverge samples for SNP calling 
# options = ["highcov", "lowcov"]
COV = ["highcov"]

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
        ## Subsample high-coverage bams
#        expand("data/final_bams/lowcov/IN{high}.subsample.dedup.bam", high = HIGH),
#        expand("reports/bamqc/final_bams/lowcov/{bam}_stats/genome_results.txt", bam = BAM),
        ## nQuire
#        hist = expand("data/nQuire/{unknown}_denoised.hist", unknown = UNKNOWN),
#        r = "data/nQuire/nquire_results_denoised.txt",
        ## Calculate AB - not effective
#        expand("data/ab/vcf/{test}.vcf", test = TEST),
#        expand("data/ab/table/{test}.AD.table", test = TEST),
#         "data/ab/table/INDE_PCRfree_Loretta_SAL-5_6_GCACGGAC_Andropogon_gerardii.AB.10k.estimate.txt"
        ## SNP Calling
#        haplocaller = expand("data/gvcf/subsample/{vcf}.p{ploidy}.subsample.g.vcf.gz", zip, ploidy = PLOIDY, vcf = VCF),
#        expand("data/gvcf/subsample/Ag.subsample.{chrom}.vcf.gz", chrom = CHROM)
#        scatter = expand("data/processed/scattered_intervals/{count}-scattered.interval_list", count = INTERVALS)
#        jointgeno = expand("data/raw/vcf_bpres/subsample/{window}.subsample.raw.vcf.gz", window = INTERVALS)
#        index = expand( "data/final_bams/highcov/{bam}.dedup.bam.bai", bam = BAM),
        mpileup = expand("data/vcf/highcov/all.AG.highcov.{chr}.raw.vcf.gz", chr = CHROM)
        ## Filtering
#        diag = expand("reports/filtering/subsample/Ag.subsample.{chrom}.table", chrom = CHOM)

# Rules
#include: "rules/mapping.smk"
#include: "rules/process_bam.smk"
#include: "rules/determine_ploidy.smk"
#include: "rules/calc_AB.smk"
include: "rules/calling.smk"
#include: "rule/filtering.smk"
