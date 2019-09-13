import config

# Different SAMPLE names
#SAMPLE = glob_wildcards("/group/jrigrp10/andropogon_shortreads/{sample}_R1.fastq.gz").sample
#SAMPLE = glob_wildcards("data/sorted_bam/{sample}.merge.bam").sample
#print(SAMPLE)

##### Test Test #####
SAMPLE = ["INCQ_PCRfree_1-McKain_252-WiDiv_43_CCGCGGTT_Andropogon_gerardii", "INCR_PCRfree_Loretta_12mi-3_25_TTATAACC_Andropogon_gerardii", "INCS_PCRfree_Loretta_CDB-2_37_GGACTTGG_Andropogon_gerardii"]
print(SAMPLE)

# Rule all describes the final output of the pipeline
rule all:
    input:
        ## Aligning reads
        #expand("data/sorted_bam/{sample}.sorted.bam", sample=SAMPLE),
        #expand("data/coverage/{sample}.cov", sample=SAMPLE),
        ## Processing the genome
        #expand("data/interm/mark_dups/{sample}.dedup.bam", sample = SAMPLE),
        #expand("data/interm/mark_dups/{sample}.dedup.bai", sample = SAMPLE),
        #expand("qc/mark_dup/{sample}_metrics.txt", sample = SAMPLE),
        ## Prepping genome
        #config.ref_dict,
        #config.ref_fai,
        ## Calling SNPs
        expand("data/gvcf/{sample}.g.vcf", sample = SAMPLE),
        config.joint_out

# Rules
#include: "rules/mapping.smk"
#include: "rules/coverage.smk"
include: "rules/process_bam.smk"
include: "rules/calling_snps.smk"
