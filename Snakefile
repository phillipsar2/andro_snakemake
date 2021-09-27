import config

# Different SAMPLE names
SAMPLE = glob_wildcards("/group/jrigrp10/andropogon_shortreads/{sample}.merge.R1.fastq.gz").sample
#SAMPLE = glob_wildcards("/group/jrigrp10/andropogon_shortreads/{sample}_1.fq.gz").sample
#SAMPLE = glob_wildcards("data/interm/mark_dups/{sample}.dedup.bam").sample
#SAMPLE = glob_wildcards("data/bcfs/{sample}.vcf").sample
#print(SAMPLE)


# BAM names encompass all samples regardless of their fastq pattern
BAM = glob_wildcards("data/interm/mark_dups/{bam}.dedup.bam").bam
#print(BAM)

## This one right below is the good shit for 11/4 work
#SAMPLE = glob_wildcards("/group/jrigrp10/andropogon_shortreads/{sample}.merge.R1.fastq.gz").sample
#SAMPLE = glob_wildcards("/group/jrigrp6/andropogon_sequence/data/bcfs/{sample}.g.vcf.gz").sample
#SAMPLE = glob_wildcards("/group/jrigrp10/andropogon_shortreads/newvcfs_fixed/{sample}.g.vcf.gz").sample

# Number of intervals for GATK SplitIntervals
#INTERVALS = ["{:04d}".format(x) for x in list(range(200))]

#SAMPLE = glob_wildcards("data/sorted_bam/{sample}.merge.bam").sample
#print(SAMPLE)

##### Test Test #####
#SAMPLE = ["INEK_PCRfree_Loretta_BB-5_22_GGTCACGA_Andropogon_gerardii", "INDZ_PCRfree_Loretta_FUL-5_35_GATCTATC_Andropogon_gerardii", "INDI_PCRfree_Loretta_SAL-6_7_TCCAACGC_Andropogon_gerardii", "INEA_PCRfree_Loretta_WAL-5_13_AGCTCGCT_Andropogon_gerardii", "INDJ_PCRfree_4-2-Oct2-plant_1_48_TACCGAGG_Andropogon_gerardii", "INEB_PCRfree_CAM807_51_ACACTAAG_Andropogon_gerardii", "INDK_PCRfree_Loretta_12mi-7_29_CGTTAGAA_Andropogon_gerardii", "INEC_PCRfree_Loretta_BB-3_20_GTGTCGGA_Andropogon_gerardii", "INDL_PCRfree_Loretta_FUL-2_32_AGCCTCAT_Andropogon_gerardii", "INED_PCRfree_Loretta_Ful-1_31_TTCCTGTT_Andropogon_gerardii", "INDM_PCRfree_Loretta_WAL-2_10_GATTCTGC_Andropogon_gerardii"] 
#print(SAMPLE)
#SAMPLE = ["INEE_PCRfree_Loretta_WAL-6_14_CCTTCACC_Andropogon_gerardii"]

# Rule all describes the final output of the pipeline
rule all:
    input:
        ## Aligning reads
        expand("data/interm/mark_dups/{sample}.dedup.bam", sample = SAMPLE),
        expand("reports/bamqc/{bam}_stats/qualimapReport.html", bam = BAM),       
        ## SNP Calling
#        expand("data/mpileup/{sample}.mpileup", sample = SAMPLE),
#        expand("data/bcfs/{sample}.vcf", sample = SAMPLE),
#        expand("data/bcfs/{sample}.vcf.gz", sample = SAMPLE),
#        expand("data/bcfs/{sample}.new.vcf.gz.csi", sample = SAMPLE),
#        expand("data/bcfs/{sample}.fixed.vcf", sample = SAMPLE),
#        expand("data/gvcf/{sample}.g.vcf.gz", sample = SAMPLE),
#        expand("/group/jrigrp10/andropogon_shortreads/newvcfs_fixed/{sample}.g.vcf.gz", sample = SAMPLE),
#        expand("/group/jrigrp10/andropogon_shortreads/newvcfs_fixed/{sample}.g.vcf.gz.tbi",sample = SAMPLE),
#        "/group/jrigrp10/andropogon_shortreads/newvcfs_fixed/combined.g.vcf.gz"
#        "data/raw/vcf/andro.raw.snps.indels.vcf",
         #config.joint_out

# Rules
include: "rules/mapping.smk"
#include: "rules/coverage.smk"
#include: "rules/process_bam.smk"
#include: "rules/calling_snps.smk"
#include: "rule/updog_genotyping.smk"
