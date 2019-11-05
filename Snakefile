import config

# Different SAMPLE names
#SAMPLE = glob_wildcards("/group/jrigrp10/andropogon_shortreads/{sample}_R1.fastq.gz").sample
SAMPLE = glob_wildcards("/group/jrigrp10/andropogon_shortreads/{sample}_L1_R1.fastq.gz").sample
#SAMPLE = glob_wildcards("data/sorted_bam/{sample}.merge.bam").sample
#print(SAMPLE)

##### Test Test #####
#SAMPLE = ["INEK_PCRfree_Loretta_BB-5_22_GGTCACGA_Andropogon_gerardii", "INDZ_PCRfree_Loretta_FUL-5_35_GATCTATC_Andropogon_gerardii", "INDI_PCRfree_Loretta_SAL-6_7_TCCAACGC_Andropogon_gerardii", "INEA_PCRfree_Loretta_WAL-5_13_AGCTCGCT_Andropogon_gerardii", "INDJ_PCRfree_4-2-Oct2-plant_1_48_TACCGAGG_Andropogon_gerardii", "INEB_PCRfree_CAM807_51_ACACTAAG_Andropogon_gerardii", "INDK_PCRfree_Loretta_12mi-7_29_CGTTAGAA_Andropogon_gerardii", "INEC_PCRfree_Loretta_BB-3_20_GTGTCGGA_Andropogon_gerardii", "INDL_PCRfree_Loretta_FUL-2_32_AGCCTCAT_Andropogon_gerardii", "INED_PCRfree_Loretta_Ful-1_31_TTCCTGTT_Andropogon_gerardii", "INDM_PCRfree_Loretta_WAL-2_10_GATTCTGC_Andropogon_gerardii"] 
#print(SAMPLE)

# Rule all describes the final output of the pipeline
rule all:
    input:
        ## Aligning reads
        #expand("data/sorted_bam/{sample}.sorted.bam", sample=SAMPLE),
        #expand("data/coverage/{sample}.cov", sample=SAMPLE),
        ## Processing the genome
        expand("data/interm/mark_dups/{sample}.dedup.bam", sample = SAMPLE),
        expand("data/interm/mark_dups/{sample}.dedup.bai", sample = SAMPLE),
        expand("qc/mark_dup/{sample}_metrics.txt", sample = SAMPLE),
        ## Prepping genome
        #config.ref_dict,
        #config.ref_fai,
        ## Calling SNPs
        #expand("data/gvcf/{sample}.g.vcf", sample = SAMPLE),
        #config.joint_out

# Rules
#include: "rules/mapping.smk"
#include: "rules/coverage.smk"
include: "rules/process_bam.smk"
#include: "rules/calling_snps.smk"
