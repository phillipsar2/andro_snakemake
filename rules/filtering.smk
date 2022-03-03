# Index vcfs
rule index_vcfs:
    input:
        "data/vcf/lowcov/all.AG.lowcov.{chr}.raw.vcf.gz"
    output:
        "data/vcf/lowcov/all.AG.lowcov.{chr}.raw.vcf.gz.tbi"
    run:
        shell("bcftools index -t {input}")

# (1) Extract SNPs from each gvcf

rule get_snps:
    input:
        ref = config.ref,
        vcf = "data/vcf/lowcov/all.AG.lowcov.{chr}.raw.vcf.gz"
    output:
         "data/raw/vcf_bpres/lowcov/all.AG.lowcov.{chr}.raw.snps.vcf.gz"
    run:
        shell("gatk SelectVariants \
        -R {input.ref} \
        -V {input.vcf} \
        -select-type SNP \
        -O {output}")


# (2) Filtering diagnostics
# Extract variant quality scores
# https://evodify.com/gatk-in-non-model-organism/

rule diagnostics:
    input:
#        vcf = "data/raw/vcf_bpres/subsample/{vcf}.p{ploidy}.subsample.raw.snps.vcf"
        vcf = "data/raw/vcf_bpres/lowcov/all.AG.lowcov.{chr}.raw.snps.vcf.gz",
        ref = config.ref
    output:
        "reports/filtering/lowcov/all.AG.lowcov.{chr}.table"
    shell:
        """
        gatk VariantsToTable \
        -R {input.ref} \
        -V {input.vcf} \
        -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR -F AD \
        -O {output}
        """

# (3) Hard filter SNPs
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531112?id=23216#2
# https://gatk.broadinstitute.org/hc/en-us/articles/360037499012?id=3225


# Hard filter for mapping quality and base quality

rule filter_snps:
    input:
        ref = config.ref,
        vcf = "data/raw/vcf_bpres/lowcov/all.AG.lowcov.{chr}.raw.snps.vcf.gz"
    output:
        "data/processed/filtered_snps_bpres/lowcov/all.AG.lowcov.{chr}.filtered.snps.vcf"
    run:
        shell("gatk VariantFiltration \
        -V {input.vcf} \
        -filter \"QUAL < 30.0\" --filter-name \"QUAL30\" \
        -filter \"MQ < 30.0\" --filter-name \"MQ30\" \
        -O {output}")


# (4) Filter SNPs to only biallelic sites and exclude the sites that failed the hard filter

rule filter_nocall:
    input:
        ref = config.ref,
        vcf = "data/processed/filtered_snps_bpres/lowcov/all.AG.lowcov.{chr}.filtered.snps.vcf"
    output:
        "data/processed/filtered_snps_bpres/lowcov/all.AG.lowcov.{chr}.filtered.nocall.vcf"
    run: 
        shell("gatk SelectVariants -V {input.vcf} --exclude-filtered true  --restrict-alleles-to BIALLELIC -O {output}")
        

# (5) Evaluate depth across samples to set DP filter

rule depth:
    input:
        vcf = "data/processed/filtered_snps_bpres/lowcov/all.AG.lowcov.{chr}.filtered.nocall.vcf",
        ref = config.ref
    output:
        dp = "reports/filtering/depth/lowcov/all.AG.lowcov.{chr}.filtered.nocall.table"
    run:
        shell("gatk VariantsToTable \
        -R {input.ref} \
        -V {input.vcf} \
        -F CHROM -F POS \
        -GF DP \
        -O {output.dp}")


# (6) Filter by genotype depth and missingness
rule filter_depth:
    input:
        vcf = "reports/filtering/depth/lowcov/all.AG.lowcov.{chr}.filtered.nocall.table",
        p = "0.99",
        miss = "0.2",
    output:
    script:
        "scripts/genoDPfilter.R"

# needs a tab deliminated list of file containing regions to select
rule keep_snps:
    input:
        vcf = "data/processed/filtered_snps_bpres/lowcov/all.AG.lowcov.{chr}.filtered.nocall.vcf",
        snps = 
    output:
        "data/processed/filtered_snps_bpres/lowcov/all.AG.lowcov.{chr}.filtered.{p}.{miss}.snps.vcf.gz"
    shell:
        """
        bgzip {input.vcf}
        tabix -p vcf {input.vcf}
        bcftools view -R {input.snps} -Oz -o {output} {input.vcf}
        """

#        shell("tabix -p vcf {output}")

# Combine interval vcfs with bcftools

#rule combine_vcfs:
#    input:
#        expand("data/processed/filtered_snps_bpres/subsample/{vcf}.p{ploidy}.subsample.filtered.3dp100.nocall.snps.vcf.gz", interval = INTERVALS, REF=REF)
#    output:
#        "data/processed/filtered_snps_bpres/{REF}/oryza.{REF}.vcf.gz"
#    run:
#        shell("bcftools concat {input} -Oz -o {output}")