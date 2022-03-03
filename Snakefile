import config
import pandas as pd

# Different SAMPLE names
#SAMPLE = glob_wildcards("/group/jrigrp10/andropogon_shortreads/{sample}.merge.R1.fastq.gz").sample
#SAMPLE = glob_wildcards("/group/jrigrp10/andropogon_shortreads/{sample}_1.fq.gz").sample
SAMPLE = glob_wildcards("/group/jrigrp10/andropogon_shortreads/ucd_seq/{sample}_R1_001.fastq.gz").sample
#print(SAMPLE)

# BAM names encompass all samples regardless of their fastq pattern
#BAM = glob_wildcards("data/interm/mark_dups/{bam}.dedup.bam").bam
BAM = glob_wildcards("data/final_bams/lowcov/{bam}.dedup.bam").bam
#BAM = glob_wildcards("data/final_bams/highcov/{bam}.dedup.bam").bam
#print(BAM)

# MERGE contains a list of the bams that belong to each genotype (GENO)  so they can be merged
#file = pd.read_csv("bams_to_merge.tsv", sep = "\t", header = 0)
file = pd.read_csv("bams_to_merge2.tsv", sep = "\t", header = 0)
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
#        expand("reports/bamqc/merged/{geno}.{merge_A}.{merge_B}_stats/genome_results.txt",zip, merge_A = MERGE_A, merge_B = MERGE_B, geno = GENO)
        ## Subsample high-coverage bams
#        expand("data/final_bams/lowcov/IN{high}.subsample.dedup.bam", high = HIGH),
#        expand("reports/bamqc/final_bams/lowcov/{bam}_stats/genome_results.txt", bam = BAM),
        ## nQuire
#        hist = expand("data/nQuire/{unknown}_denoised.hist", unknown = UNKNOWN),
#        r = "data/nQuire/nquire_results_denoised.txt",
        ## SNP Calling
#        index = expand("data/final_bams/lowcov/{bam}.dedup.bam.bai", bam = BAM),
#        mpileup = expand("data/vcf/highcov/all.AG.highcov.{chr}.raw.vcf.gz", chr = CHROM)
#        mpileup = expand("data/vcf/lowcov/ucd.lowcov.{chr}.raw.vcf.gz", chr = CHROM)
        ## Filtering
#        index = expand("data/vcf/lowcov/all.AG.lowcov.{chr}.raw.vcf.gz.tbi", chr = CHROM),
#        snp = expand("data/raw/vcf_bpres/lowcov/all.AG.lowcov.{chr}.raw.snps.vcf.gz", chr = CHROM),
#        diag = expand("reports/filtering/lowcov/all.AG.lowcov.{chr}.table", chr = CHROM)
#        hf = expand("data/processed/filtered_snps_bpres/lowcov/all.AG.lowcov.{chr}.filtered.nocall.vcf", chr = CHROM)
        dp_diag = expand("reports/filtering/depth/lowcov/all.AG.lowcov.{chr}.filtered.nocall.table", chr = CHROM)

# Rules
#include: "rules/mapping.smk"
#include: "rules/process_bam.smk"
#include: "rules/determine_ploidy.smk"
#include: "rules/calling.smk"
include: "rules/filtering.smk"
