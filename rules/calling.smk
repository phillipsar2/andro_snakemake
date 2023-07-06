# index

rule index:
    input:
        "data/final_bams/lowcov/{bam}.dedup.bam"
#        "data/final_bams/highcov/{bam}.dedup.bam"
    output:
        "data/final_bams/lowcov/{bam}.dedup.bam.bai"
#        "data/final_bams/highcov/{bam}.dedup.bam.bai"
    shell:
        """
        samtools index {input}
        """

## Potentially calling snps with bcftools because we just want to get read counts
# calling SNPs for the population panel then calculating quality and depth statistics for filtering
rule mpileup:
    input:
        ref = config.ref,
#        bamlist = "data/final_bams/lowcov/lowcov.bamlist.txt"
        bamlist = "data/final_bams/highcov/highcov.bamlist.txt"
#        bamlist = "data/final_bams/lowcov/ucd.bamlist.txt"
    output:
#        vcf = "data/vcf/lowcov/all.AG.lowcov.{chr}.raw.vcf.gz"
        vcf = "data/vcf/highcov/all.AG.highcov.{chr}.raw.vcf.gz"
#        vcf = "data/vcf/lowcov/ucd.lowcov.{chr}.raw.vcf.gz"
    params:
        chr = "{chr}"
    threads:
        8
    run:
       # default only sites with max 250 reads considered at each positin, this is way above the max coverage
       # -v option asks to output variant sites only (this is sufficient for the analyses we want to run)
       # -r output for only the given region
        shell("bcftools mpileup -Ou -f {input.ref} -b {input.bamlist} -r {params.chr} \
        --annotate FORMAT/AD,FORMAT/DP --threads {threads} | \
        bcftools call -mv -Oz -o {output.vcf}")
        shell("bcftools index -t {output}")

# select samples to keep from the vcfs
rule drop_samples:
    input:
        samp_list = "data/vcf/lowcov/all.AG.samples",
        vcf = "data/vcf/lowcov/all.AG.lowcov.{chr}.raw.vcf.gz"
    output:
        "data/vcf/lowcov/AG.temp.lowcov.{chr}.vcf.gz"
    shell:
        """
        bcftools view -S {input.samp_list} -Oz -o {output} {input.vcf}
        bcftools index -t {output}
        """


# merge vcfs that samples were dropped from and combine with new ucd vcfs
rule merge_vcfs:
    input:
        ucd = "data/vcf/lowcov/ucd.lowcov.{chr}.raw.vcf.gz",
        ag = "data/vcf/lowcov/AG.temp.lowcov.{chr}.vcf.gz"
    output:
        "data/vcf/lowcov/all.AG.lowcov.{chr}.raw.vcf.gz"
    shell:
        """
        bcftools merge {input.ucd} {input.ag} -Oz -o {output}
        """
