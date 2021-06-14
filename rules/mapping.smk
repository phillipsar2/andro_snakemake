# Align reads to the reference genome
rule bwa_map:
    input:
        ref = config.ref,
        r1 = "/group/jrigrp10/andropogon_shortreads/{sample}_1.fq.gz",
        r2 = "/group/jrigrp10/andropogon_shortreads/{sample}_2.fq.gz",
#        r1 = "/group/jrigrp10/andropogon_shortreads/{sample}.merge.R1.fastq.gz",
#        r2 = "/group/jrigrp10/andropogon_shortreads/{sample}.merge.R2.fastq.gz"
#
    output:
        temp("data/interm/mapped_bam/{sample}.mapped.bam"),
#    log:
#        "logs/bwa_mem/{sample}.log",
    shell:
        "bwa-mem2 mem -t 8 {input.ref} {input.r1} {input.r2} |"
        "samtools view -Sb > {output}"

# Takes the input file and stores a sorted version in a different directory.
rule samtools_sort:
    input:
        "data/interm/mapped_bam/{sample}.mapped.bam"
    output:
        temp("data/interm/sorted_bam/{sample}.sorted.bam"),
    params:
        tmp = "/scratch/aphillip/sort_bam/{sample}"
    run:
        shell("mkdir -p {params.tmp}")
        shell("samtools sort -T {params.tmp} {input} > {output}")
        shell("rm -rf {params.tmp}")

rule add_rg:
    input:
        "data/interm/sorted_bam/{sample}.sorted.bam"
    output:
        bam = temp(touch("data/interm/addrg/{sample}.rg.bam"))
    params:
        tmp = "/scratch/aphillip/addrg/{sample}",
        sample = "{sample}"
    run:
        shell("mkdir -p {params.tmp}")
        shell("gatk AddOrReplaceReadGroups \
        -I={input} \
        -O={output.bam} \
        -RGID=4 \
        -RGLB=lib1 \
        -RGPL=illumina \
        -RGPU=unit1 \
        -RGSM={params.sample} \
        --TMP_DIR {params.tmp} \
        --CREATE_INDEX=true")
        shell("rm -rf {params.tmp}")

rule mark_dups:
    input:
        "data/interm/addrg/{sample}.rg.bam"
    output:
        bam = "data/interm/mark_dups/{sample}.dedup.bam",
        metrics = "qc/mark_dup/{sample}_metrics.txt"
    params:
        tmp = "/scratch/aphillip/mark_dups/{sample}"
    run:
        # Create a scratch directory
        shell("mkdir -p {params.tmp}")
        # Input bam file to output marked records. Assume bam file has been sorted. Direct to a temporary storage file (scratch).
        shell("gatk MarkDuplicates \
        -I={input} \
        -O={output.bam} \
        --METRICS_FILE={output.metrics} \
        --CREATE_INDEX=true \
        -MAX_FILE_HANDLES=1000 \
        --ASSUME_SORT_ORDER=coordinate \
        --TMP_DIR={params.tmp}")
        # Remove scratch directory
        shell("rm -rf {params.tmp}")

# Quality metrics with qualimap
rule bamqc:
    input:
        "data/interm/mark_dups/{sample}.dedup.bam"
    output:
        "reports/bamqc/{sample}_stats/qualimapReport.html"
    params:
        dir = "reports/bamqc/{sample}_stats"
    run: 
        shell("qualimap bamqc \
        -bam {input} \
        -nt 8 \
        -nr 100000 \
        -outdir {params.dir} \
        -outformat HTML \
        --skip-duplicated \
        --java-mem-size=64G")

### No longer need to realign INDELs as HaplotypeCaller takes care of it ###
# https://gatkforums.broadinstitute.org/gatk/discussion/11455/realignertargetcreator-and-indelrealigner
        
