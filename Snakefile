# Mapping Andropogon reads to the reference genomes
SAMPLES = glob_wildcards("/home/aphillip/Andropogoneae/Andro_sequences/Andropogonreads/{sample}_R1.fastq.gz").sample

# Rule all describes the final output of the pipeline
rule all:
    input:
        # sorted bams - sorted bam files for each sample
        expand("/home/aphillip/Andro_snakemake/data/sorted_reads/{sample}.sorted.bam", sample=SAMPLES),
#        expand("/home/aphillip/Andro_snakemake/data/coverage/{sample}.cov", sample=SAMPLES)

# Align reads to the reference genome
rule bwa_map:
    input:
        ref = "/home/aphillip/Andropogoneae/Andro_sequences/ANDRO.asmAll.contigs.fasta.gz",
        r1 = "/home/aphillip/Andropogoneae/Andro_sequences/Andropogonreads/{sample}_R1.fastq.gz",
        r2 = "/home/aphillip/Andropogoneae/Andro_sequences/Andropogonreads/{sample}_R2.fastq.gz",
    output:
        temp("/home/aphillip/Andro_snakemake/data/interm/mapped_reads/{sample}.mapped.bam"),
    log:
        "/home/aphillip/Andro_snakemake/logs/bwa_mem/{sample}.log",
    #threads: 24
    shell:
#        shell("bwa index {input.ref}")
        "(bwa mem -t 24 {input.ref} {input.r1} {input.r2} |"
        "samtools view -Sb > {output}) 2> {log}"

# Takes the input file and stores a sorted version in a different directory.
rule samtools_sort:
    input:
        "/home/aphillip/Andro_snakemake/data/interm/mapped_reads/{sample}.mapped.bam",
    output:
        "/home/aphillip/Andro_snakemake/data/sorted_reads/{sample}.sorted.bam",
    shell:
        "samtools sort {input} > {output}"

# Calculate average coverage across each sorted bam file
#rule find_coverage:
#    input:
#        "/home/aphillip/Andro_snakemake/data/sorted_reads/{sample}.sorted.bam",
#    output:
#        "/home/aphillip/Andro_snakemake/data/coverage/{sample}.cov",
#    shell:
#        # Have to double up on brackets, backslashes, apostrophes, etc.
#        "samtools depth {input} | "
#        "awk '{{sum+=$3;sumsq+=$3*$3}} END {{ print ""Average ="",sum/5105873300; print ""Stdev = "",sqrt(sumsq/5105873300 - (sum/5105873300)**2)}}' > {output} "

