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
        "data/ebg/lowcov/total_reads.{chr}.{ploidy}.txt",
        "data/ebg/lowcov/alt_reads.{chr}.{ploidy}.txt",
        "data/ebg/lowcov/error_rate.{chr}.{ploidy}.txt",
        "data/ebg/lowcov/samples.{chr}.{ploidy}.txt",
        "data/ebg/lowcov/snp_positions.{chr}.{ploidy}.txt"
    params:
        chr = "{chr}",
        ploidy = "{ploidy}"
    shell:
        "Rscript scripts/vcf2ADmatrix.R {input} --chr {params.chr} --ploidy {params.ploidy}"

# (3) Run EBG

rule ebg:
    input:
        tot = "data/ebg/lowcov/total_reads.{chr}.{ploidy}.txt",
        alt = "data/ebg/lowcov/alt_reads.{chr}.{ploidy}.txt",
        e = "data/ebg/lowcov/error_rate.{chr}.{ploidy}.txt",
        samp = "data/ebg/lowcov/samples.{chr}.{ploidy}.txt",
        loci = "data/ebg/lowcov/snp_positions.{chr}.{ploidy}.txt"
    output:
        "data/ebg/lowcov/{chr}.{ploidy}-PL.txt"
    params:
        p = "2",
        pre = "data/ebg/lowcov/{chr}.{ploidy}"
    shell:
        """
        samples=`wc -l {input.samp}`
        loci=`wc -l {input.loci}`
        ~/toolsfordayz/polyploid-genotyping/ebg/ebg diseq \
        -t {input.tot} \
        -a {input.alt} \
        -e {input.e} \
        -p {params.p} \
        -n $samples \
        -l $loci \
        --iters 1000 \
        --prefix {params.pre}
        """
