# Mapping Andropogon reads to the reference genomes
import config

# Different SAMPLE names
#SAMPLE = glob_wildcards("/group/jrigrp10/andropogon_shortreads/{sample}_R1.fastq.gz").sample
SAMPLE = glob_wildcards("data/sorted_bam/{sample}.merge.bam").sample
#print(SAMPLE)

# Rule all describes the final output of the pipeline
rule all:
    input:
        # sorted bams - sorted bam files for each sample
     #expand("data/sorted_bam/{sample}.sorted.bam", sample=SAMPLE),
     expand("data/coverage/{sample}.cov", sample=SAMPLE)

# Rules
#include: "rules/mapping.smk"
include: "rules/coverage.smk"
#include: "rules/calling_snps.smk"
