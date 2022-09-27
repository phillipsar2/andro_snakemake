# Call genotype likelihoods with EBG

# (1) split vcf by ploidy

rule split_ploidy:
    input:
        vcf = "data/processed/filtered_snps_bpres/lowcov/all.AG.lowcov.{chr}.filtered.99.20.snps.vcf.gz",
        geno = "{ploidy}.lowcov.samples.csv"
    output:
        temp("data/processed/filtered_snps_bpres/lowcov/AG.lowcov.{chr}.{ploidy}.snps.vcf")
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
# need to manually change p from 2 to 3

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

# (4) convert EBG output to matrix format

rule gl_mat:
    input:
        tot = "data/ebg/lowcov/total_reads.{chr}.6x.txt",
        tot2 = "data/ebg/lowcov/total_reads.{chr}.9x.txt",
        geno = "data/ebg/lowcov/samples.{chr}.6x.txt",
        geno2 = "data/ebg/lowcov/samples.{chr}.9x.txt",
        snps = "data/ebg/lowcov/snp_positions.{chr}.6x.txt",
        pl_hex = "data/ebg/lowcov/{chr}.6x-PL.txt",
        pl_enn = "data/ebg/lowcov/{chr}.9x-PL.txt"
    output:
        "data/ebg/lowcov/{chr}-GL.txt"
    shell:
        """
        Rscript scripts/ebg2glmatrix_2ploidy.R --tot {input.tot} -p 2 --tot2 {input.tot2} --ploidy2 3 --geno {input.geno} --geno2 {input.geno2} -s {input.snps} {input.pl_hex} {input.pl_enn}
        """

# (5a) Merge GL.txt files and ranomly select 10k snps
# double-check no header lines were randomly selected afterwards

#rule 10k_snps:
#    input:
#        expand("data/ebg/lowcov/{chr}-GL.txt", chr = CHROM)
#    output:
#        "data/ebg/lowcov/10k_lowcov-GL.txt"
#    shell:
#        """
#        head -n1 data/ebg/lowcov/Chr01A-GL.txt > {output}
#        shuf -n 10000 <(cat {input}) >> {output}
#        """

# for high cov:
# head -n1 data/processed/filtered_snps_bpres/highcov/all.AG.highcov.Chr01A.PL.txt > data/processed/filtered_snps_bpres/highcov/10k_highcov-PL.txt
# shuf -n 10000 <(cat data/processed/filtered_snps_bpres/highcov/*.PL.txt) >> data/processed/filtered_snps_bpres/highcov/10k_highcov-PL.txt

# (5b) Select sites with less than X% missing data
rule filt_miss:
    input:
        "data/ebg/lowcov/genoliks/{chr}-GL.txt"
    output:
#        "data/ebg/lowcov/{chr}.miss5-GL.txt"
        "data/ebg/lowcov/genoliks/{chr}.miss20-GL.txt"
    shell:
        "Rscript scripts/strict_missing_data.R {input}"


# ----------- High coverage -------------

# (1) grab GL from the vcfs
rule grab_gl:
    input:
        vcf = "data/processed/filtered_snps_bpres/{cov}/all.AG.{cov}.{chr}.filtered.99.20.snps.vcf.gz",
        ref = config.ref
    output:
        "data/processed/filtered_snps_bpres/{cov}/all.AG.{cov}.{chr}.PL.txt"
    shell:
        """
        gatk IndexFeatureFile -I {input.vcf}
        gatk VariantsToTable \
        -R {input.ref} \
        -V {input.vcf} \
        -F CHROM -F POS -GF PL \
        -O {output}
        """

