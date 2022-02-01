# index

rule index:
    input:
#        "data/final_bams/lowcov/{bam}.dedup.bam"
        "data/final_bams/highcov/{bam}.dedup.bam"
    output:
#        "data/final_bams/lowcov/{bam}.dedup.bam.bai"
        "data/final_bams/highcov/{bam}.dedup.bam.bai"
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
    output:
#        vcf = "data/vcf/lowcov/all.AG.lowcov.{chr}.raw.vcf.gz"
        vcf = "data/vcf/highcov/all.AG.highcov.{chr}.raw.vcf.gz"
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




# Calling variant sites using GATK
# Follows GATK best practices

# (1) Call SNPs and genotype likelihoods with HaplotypeCaller
# keep invariant sites http://seqanswers.com/forums/showthread.php?t=43006
## include --emitRefConfidence (-ERC) GVCF or BP_RESOLUTION in HaplotypeCaller to get genotype likelihoods
rule haplotype_caller:
    input:
        ref = config.ref, 
#        bam = "data/interm/mark_dups/{vcf}.dedup.bam"
        bam = "data/final_bams/lowcov/{vcf}.subsample.dedup.bam"
    output:
#        outdir = "data/gvcf/{vcf}.p{ploidy}.g.vcf.gz"
        outdir = "data/gvcf/subsample/{vcf}.p{ploidy}.subsample.g.vcf.gz"
    params:
        regions = config.contig_list,
        ploidy = "{ploidy}"
    run:
        shell("samtools index {input.bam}")
        shell("gatk --java-options \"-Xms96G -Xmx96G\" \
        HaplotypeCaller \
        --input {input.bam} \
        --output {output.outdir} \
        --reference {input.ref} \
        --G StandardAnnotation \
        -G AS_StandardAnnotation \
        -L {params.regions} \
        -ploidy {params.ploidy} \
        -ERC BP_RESOLUTION \
        --native-pair-hmm-threads 8")


### Merge vcfs horizontally, seperate by chromosome, and grab SNPs & called genotypes only
rule merge_vcfs:
    input:
        expand("data/gvcf/subsample/{vcf}.p{ploidy}.subsample.g.vcf.gz", zip, vcf = VCF, ploidy = PLOIDY)
    output:
        "data/gvcf/subsample/Ag.subsample.{chrom}.vcf.gz"
    params:
        chrom = "{chrom}"
    threads:
        8
    shell:
        """
        bcftools merge --file-list {input} -R {params.chrom} --threads {threads} -Ou | \
        bcftools view -v snps -U -Oz -o {output}
        """



# ----- not using below


# Scatter reference into intervals using SplitIntervals. Set intervals in Snakefile.
# https://gatk.broadinstitute.org/hc/en-us/articles/360036348372-SplitIntervals
# ``--scatter-count n` splits reference into n intervals
rule split_intervals:
    input:
        ref = config.ref
    output:
        int = expand("data/processed/scattered_intervals/{count}-scattered.interval_list", count = INTERVALS)
    params:
        regions = config.contig_list,
        dir = "data/processed/scattered_intervals/"
    run:
        shell("gatk SplitIntervals -R {input.ref} -L {params.regions} --scatter-count 200 -O {params.dir}")

# bgzip vcfs and tbi index
rule bgzip:
    input:
        "data/gvcf/subsample/{vcf}.p{ploidy}.subsample.g.vcf.gz"
    output:
        "data/gvcf/subsample/{vcf}.p{ploidy}.subsample.g.vcf.gz.tbi"
    shell:
        """
        gunzip -c {input} | bgzip > {input}
        bcftools index -t {input}
        """


# Combine GVCFs with GenomicsDBImport
# sample_map file is samplename with path/to/gvcf in the following line
# scattered interval file created by SpitIntevals

rule combine_gvcfs:
    input:
#        gvcfs = expand("data/gvcf/{vcf}.p{ploidy}.g.vcf.gz", zip, vcf = VCF, ploidy = PLOIDY),
        gvcfs = expand("data/gvcf/subsample/{vcf}.p{ploidy}.subsample.g.vcf.gz", zip, vcf = VCF, ploidy = PLOIDY),
        region = "data/processed/scattered_intervals/{window}-scattered.interval_list",
        map = config.sample_map
    output:
#        directory("data/interm/combined_database_bpres/{count}")
        directory("data/interm/combined_database_bpres/subsample/{window}")
    params:
        tmp = "/scratch/aphillip/genomicsdbimport/{window}"
    shell:
        """
        mkdir -p {params.tmp}
        gatk --java-options \"-Xmx90g -Xms90g\" \
        GenomicsDBImport \
        --genomicsdb-workspace-path {output} \
        --batch-size 50 \
        --reader-threads 8 \
        --sample-name-map {input.map} \
        --intervals {input.region} --tmp-dir {params.tmp}
        rm -rf {params.tmp}
        """


rule joint_geno:
    input:
#        dir = directory("data/interm/combined_database_bpres/{interval}"),
        dir = directory(expand("data/interm/combined_database_bpres/subsample/{window}", window = INTERVALS)),
        ref = config.ref
    output:
#        "data/raw/vcf_bpres/{interval}.raw.vcf.gz"
        "data/raw/vcf_bpres/subsample/{interval}.subsample.raw.vcf.gz"
    params:
        db = "gendb://data/interm/combined_database_bpres/{window}",
        region = "data/processed/scattered_intervals/{window}-scattered.interval_list",
    run:
        shell("gatk  --java-options \"-Xms32G -Xmx32G\" \
        GenotypeGVCFs \
        -R {input.ref} \
        -V {params.db} \
        -L {params.region} \
        -new-qual \
        -G StandardAnnotation \
        -G AS_StandardAnnotation \
        --include-non-variant-sites \
        -O {output}")
