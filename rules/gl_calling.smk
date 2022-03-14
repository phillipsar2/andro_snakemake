# Call genotype likelihoods with EBG

# (1) split vcf by ploidy

rule split_ploidy:
    input:
        vcf = "data/processed/filtered_snps_bpres/lowcov/all.AG.lowcov.{chr}.filtered.99.20.snps.vcf.gz",
        geno = "{ploidy}.lowcov.samples.csv"
    output:
        "data/processed/filtered_snps_bpres/lowcov/AG.lowcov.{chr}.{ploidy}.snps.vcf"
    shell:
        """
        bcftools view -S {input.geno} --force-samples {input.vcf} > {output}
        """

# (2) Extract allele depth matrices from the vcfs and generate other EBG input files

rule vcf2ADmatrix:
    input:
        "data/processed/filtered_snps_bpres/lowcov/AG.lowcov.{chr}.{ploidy}.snps.vcf"
    output:
        "data/ebg/total_reads.{chr}.{ploidy}.txt",
        "data/ebg/alt_reads.{chr}.{ploidy}.txt",
        "data/ebg/error_rate.{chr}.{ploidy}.txt",
        "data/ebg/samples.{chr}.{ploidy}.txt",
        "data/ebg/snp_positions.{chr}.{ploidy}.txt"
    script:
        "/group/jrigrp6/andropogon_sequence/scripts/vcf2ADmatrix.R"
#    shell:
#        "Rscript scripts/genoDPfilter.R {input.vcf} -q {params.p} -m {params.miss}"
