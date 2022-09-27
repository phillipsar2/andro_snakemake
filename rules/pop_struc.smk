# Evaluating population structure with ENTROPY

# (1a) Generate ENTROPY input files from gl matrix
rule ent_in:
    input:
        "data/ebg/lowcov/10k.lowcov.miss5-GL.txt"
    output:
        "data/entropy/lowcov/10k_lowcov.miss5.qk2inds.txt",
        "data/entropy/lowcov/10k_lowcov.miss5_pntest_meangl.txt",
        "data/entropy/lowcov/10k_lowcov.miss5.mpgl"
    params:
        pre = "data/entropy/lowcov/10k_lowcov.miss5"
    shell:
        """
#        module load spack/R/4.1.1
        Rscript scripts/glmat2mpgl.R --gl {input} --out {params.pre}
        """

# (1b) Generate ENTROPY input files from PLs extracted from vcf

rule ent_in_tab:
    input:
        "data/processed/filtered_snps_bpres/highcov/10k_highcov-PL.txt"
    output:
        "data/entropy/highcov/10k_highcov.qk2inds.txt",
        "data/entropy/highcov/10k_highcov_pntest_meangl.txt",
        "data/entropy/highov/10k_highcov.mpgl"
    params:
        pre = "data/entropy/highcov/10k_highcov"
    shell:
        """
        module load spack/R/4.1.1
        Rscript scripts/gltable2mpgl.R --pl {input} --out {params.pre}
        """

# (2) Run ENTROPY
# not sure if I can actually activate a conda env in a snakemake rule?? might have to run it sepearately outside snakemake
# run with default burn-in of 1000 (-b 1000) and 10000 MCMC steps (-l 10000); store every 5 step after burn-in (-t 5)
# recommends running three times
#rule entropy:
#    input:
#        mpgl = "data/entropy/{cov}/10k_{cov}.mpgl",
#        ploidy = "data/entropy/{cov}/10k_{cov}.ploidy_inds.txt",
#        init_admix = "data/entropy/{cov}/10k_{cov}.qk{k}inds.txt",
#    output:
#        "data/entropy/{cov}/10k_{cov}.k{k}.c{chain}.hdf5"
#    params:
#        k = "{k}"
#    shell:
#        """
#        module load bio3 
#        source activate entropy-2.0
#        entropy -i {input.mpgl} \
#        -n {input.ploidy} \
#        -k {params.k} \
#        -q {input.init_admix} \
#        -t 5 \ 
#        -o {output}
#        """

# (3) Extract ENTROPY output files

rule ent_out:
    input:
        "data/entropy/{cov}/10k_{cov}.k{k}.c{chain}.hdf5"
    output:
        "data/entropy/{cov}/qest.k{k}.c{chain}.txt"
    shell:
        """
#        module load spack/R/4.1.1
        Rscript scripts/assessconvergence.R {input}
        estpost.entropy -p q -s 0 {input} -o {output}
        """


### ANGSD analyses on 6x high coverage individuals
### PCA --------

# (1) generate GL in beagle format
# doMajorMinor 4: use refence allele as major
rule angsd_beagle:
    input:
        ref = config.ref,
        bamlist = "data/final_bams/highcov/highcov.bamlist.txt"
    output:
        "data/angsd/pca/highcov.{chrom}.8dp70.beagle.gz"
    params:
        prefix = "data/angsd/pca/highcov.{chrom}.8dp70",
        chrom = "{chrom}"
    run:
        shell("angsd -GL 1 -P 15 \
        -doGlf 2 \
        -doMajorMinor 4 \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
        -minMapQ 30 -minQ 30 \
        -doCounts 1 \
        -setMinDepthInd 8 -setMaxDepthInd 70 \
        -minInd 8 \
        -ref {input.ref} -anc {input.ref} \
        -r {params.chrom} \
        -bam {input.bamlist} \
        -out {params.prefix}")

# (2) merge beagle files together
#        shell("cat <(zcat data/angsd/pca/*beagle.gz | head -n1) <(zcat data/angsd/pca/*beagle.gz | \
#        grep -v marker) > data/angsd/pca/highcov.merged.8dp70.beagle
#        shell("gzip data/angsd/pca/highcov.merged.8dp70.beagle")

# (3) grab 30k random snps
#        zcat highcov.merged.8dp70.beagle.gz | head -n1 > highcov.merged.8dp70.30k.beagle
#        shuf -n 30000 <(zcat highcov.merged.8dp70.beagle.gz) >> highcov.merged.8dp70.30k.beagle
#        gzip highcov.merged.8dp70.30k.beagle

# (4) run PCA
rule angsd_pca:
    input:
        beagle = expand("data/angsd/pca/highcov.merged.8dp70.30k.beagle.gz", chrom = CHROM)
    output:
        pca = "data/angsd/pca/highcov.merged.8dp70.30k..cov"
    run:
#       shell("module load bio3")
#       shell("source activate pcangsd")
        shell("pcangsd -beagle {input.beagle} -o data/angsd/pca/highcov.merged.8dp70.30k")

# (5) Plot in R

### PCA on low coverage dataset

# (1) grab 10k random SNPs
# grab every 9th row
# awk 'NR % 9 == 0' all.lowcov.miss5-GL.txt >> 10k.lowcov.miss5-GL.txt

# (2) Calculate allele frequencies and run PCA (see R script)
