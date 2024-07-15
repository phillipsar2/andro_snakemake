# (12) Index gvcfs
rule index_vcfs:
    input:
        "data/vcf/{cov}/all.AG.{cov}.{chr}.raw.vcf.gz"
    output:
        "data/vcf/{cov}/all.AG.{cov}.{chr}.raw.vcf.gz.tbi"
    run:
        shell("bcftools index -t {input}")

# (13) Extract SNPs from each gvcf
rule get_snps:
    input:
        ref = config.ref,
        vcf = "data/vcf/{cov}/all.AG.{cov}.{chr}.raw.vcf.gz"
    output:
         "data/raw/vcf_bpres/{cov}/all.AG.{cov}.{chr}.raw.snps.vcf.gz"
    run:
        shell("gatk SelectVariants \
        -R {input.ref} \
        -V {input.vcf} \
        -select-type SNP \
        -O {output}")


# (14) Filtering diagnostics - extract variant quality scores
# Roughly following suggestions in https://evodify.com/gatk-in-non-model-organism/
rule diagnostics:
    input:
#        vcf = "data/raw/vcf_bpres/subsample/{vcf}.p{ploidy}.subsample.raw.snps.vcf"
        vcf = "data/raw/vcf_bpres/{cov}/all.AG.{cov}.{chr}.raw.snps.vcf.gz",
        ref = config.ref
    output:
        "reports/filtering/{cov}/all.AG.{cov}.{chr}.table"
    shell:
        """
        gatk VariantsToTable \
        -R {input.ref} \
        -V {input.vcf} \
        -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR -F AD \
        -O {output}
        """

# (15) Hard filter SNPs
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531112?id=23216#2
# https://gatk.broadinstitute.org/hc/en-us/articles/360037499012?id=3225

# Hard filter for mapping quality and base quality
rule filter_snps:
    input:
        ref = config.ref,
        vcf = "data/raw/vcf_bpres/{cov}/all.AG.{cov}.{chr}.raw.snps.vcf.gz"
    output:
        "data/processed/filtered_snps_bpres/{cov}/all.AG.{cov}.{chr}.filtered.snps.vcf"
    run:
        shell("gatk VariantFiltration \
        -V {input.vcf} \
        -filter \"QUAL < 30.0\" --filter-name \"QUAL30\" \
        -filter \"MQ < 30.0\" --filter-name \"MQ30\" \
        -O {output}")


# (16) Filter SNPs to only biallelic sites and exclude the sites that failed the hard filter
rule filter_nocall:
    input:
        ref = config.ref,
        vcf = "data/processed/filtered_snps_bpres/{cov}/all.AG.{cov}.{chr}.filtered.snps.vcf"
    output:
        "data/processed/filtered_snps_bpres/{cov}/all.AG.{cov}.{chr}.filtered.nocall.vcf"
    run: 
        shell("gatk SelectVariants -V {input.vcf} --exclude-filtered true  --restrict-alleles-to BIALLELIC -O {output}")
        

# (17) Exclude clones
# Clones were identified by looking at the kinship matrix using all genotypes.
rule exclude_clones:
    input:
        vcf = "data/processed/filtered_snps_bpres/{cov}/all.AG.{cov}.{chr}.filtered.nocall.vcf",
        ref = config.ref,
        clones = "clones.list"
    output:
        vcf = "data/processed/filtered_snps_bpres/{cov}/all.AG.noclones.{cov}.{chr}.filtered.nocall.vcf",
    shell:
        """
        vcftools --vcf {input.vcf} --remove {input.clones} --stdout > {output.vcf}
        """ 

# (18) Extract genotype depth across samples to determine DP cutoff
rule depth:
    input:
#        vcf = "data/processed/filtered_snps_bpres/{cov}/all.AG.{cov}.{chr}.filtered.nocall.vcf",
        vcf = "data/processed/filtered_snps_bpres/{cov}/all.AG.noclones.{cov}.{chr}.filtered.nocall.vcf",
        ref = config.ref
    output:
#        dp = "reports/filtering/depth/{cov}/all.AG.{cov}.{chr}.filtered.nocall.table"
        dp = "reports/filtering/depth/{cov}/all.AG.noclones.{cov}.{chr}.filtered.nocall.table"
    run:
        shell("gatk VariantsToTable \
        -R {input.ref} \
        -V {input.vcf} \
        -F CHROM -F POS \
        -GF DP \
        -O {output.dp}")

# (19) Filter by genotype depth and missingness
# p = probability of the given read depth
# miss = percent missing data threshold at a site
# min = minimum depth for a genotype at a site
rule filter_depth:
    input:
#        vcf = "reports/filtering/depth/{cov}/all.AG.{cov}.{chr}.filtered.nocall.table"
        vcf = "reports/filtering/depth/{cov}/all.AG.noclones.{cov}.{chr}.filtered.nocall.table"
    output:
        ## Output used to select sites in ANGSD
#        "reports/filtering/depth/{cov}/all.AG.{cov}.{chr}.filtered.nocall.0.99_0.2.txt"
        "reports/filtering/depth/{cov}/all.AG.noclones.{cov}.{chr}.filtered.nocall.0.99_0.2.txt"
    params:
#        p = "{p}",
#        miss = "{miss}"
        p = "0.99",
        miss = "0.2",
        min = "1"
    shell:
        "Rscript scripts/genoDPfilter.R {input.vcf} -q {params.p} -m {params.miss} --min {params.min}"

