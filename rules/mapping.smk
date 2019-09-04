# Align reads to the reference genome
#SAMPLE = glob_wildcards("/group/jrigrp10/andropogon_shortreads/{sample}_R1.fastq.gz").sample
#print(SAMPLE)

rule bwa_map:
    input:
        #config.ref,
        ref = "/group/jrigrp6/andropogon_sequence/data/genome/ANDRO_contigs1278.fasta",
        r1 = "/group/jrigrp10/andropogon_shortreads/{sample}_R1.fastq.gz",
        r2 = "/group/jrigrp10/andropogon_shortreads/{sample}_R2.fastq.gz",
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
        "data/sorted_bam/{sample}.sorted.bam",
    shell:
        "samtools sort {input} > {output}"
