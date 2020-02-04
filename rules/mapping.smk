# Align reads to the reference genome
#SAMPLE = glob_wildcards("/group/jrigrp10/andropogon_shortreads/{sample}_R1.fastq.gz").sample
#print(SAMPLE)

rule bwa_map:
    input:
        #config.ref,
        ref = "/group/jrigrp6/andropogon_sequence/data/genome/newgenome/andropogon-geardii_Canu-BioNano-Hybrid_all-contigs.fasta",
        r1 = "/group/jrigrp10/andropogon_shortreads/{sample}_1.fq.gz",
        r2 = "/group/jrigrp10/andropogon_shortreads/{sample}_2.fq.gz",
#        r1 = "/group/jrigrp10/andropogon_shortreads/{sample}.merge.R1.fastq.gz",
#        r2 = "/group/jrigrp10/andropogon_shortreads/{sample}.merge.R2.fastq.gz"
    output:
        temp("data/interm/mapped_bam/{sample}.mapped.bam"),
    log:
        "logs/bwa_mem/{sample}.log",
    shell:
        "(bwa mem -t 8 {input.ref} {input.r1} {input.r2} |"
        "samtools view -Sb > {output}) 2> {log}"

# Takes the input file and stores a sorted version in a different directory.
rule samtools_sort:
    input:
        "data/interm/mapped_bam/{sample}.mapped.bam",
    output:
        temp("data/sorted_bam/{sample}.sorted.bam"),
    params:
        tmp = "/scratch/aphillip/sort_bam/{sample}"
    run:
        shell("mkdir -p {params.tmp}")
        shell("samtools sort -T {params.tmp} {input} > {output}")
        shell("rm -rf {params.tmp}")

#    output:
#        temp("data/mergensort/{sample}.merge.sorted.bam")
#    params:
#        "-m 16G"
#        tmp = "/scratch/aphillip/sort_bam/{sample}"
#    run:
#        shell("samtools sort {input} > {output}")
        #shell("rm {sample}.merge.bam")
#    wrapper:
#        "0.38.0/bio/samtools/sort"
