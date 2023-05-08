## Rules for determining if aneuploidy occurs in our samples


# (1) get windows across the genome
#  bedtools makewindows -g {config.chrom_lengths} -w 1000000 > {output.bed}

# (2) Convert bams to bed files

rule bam2bed:
    input:
        bam = "data/final_bams/lowcov/{bam}.dedup.bam"
    output:
        bed = temp("data/bedtools/bed/lowcov/{bam}.bed")
    shell:
        """
        bedtools bamtobed -i {input.bam} > {output.bed}
        """
    

# (2) Extract depth across genome in 1Mb windows for each bam file
# coverage is calculate for -a 
rule get_cov:
    input:
        ref = config.ref,
        bed = "data/bedtools/bed/lowcov/{bam}.bed",
        windows = config.windows
    output:
        cov = "data/bedtools/coverage/{bam}.1Mb.cov.txt"
    shell:
        """
        bedtools coverage -b {input.bed} -a {input.windows} > {output.cov}
        """



