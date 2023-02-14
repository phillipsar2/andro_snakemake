import random

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

###
### PCA on low coverage dataset - single read sampling approach
###

# (1) Create the sites files
rule create_sites:
    input:
        sites = "data/ebg/lowcov/genoliks/lowcov.merged.miss20.20k.positions",
        vcf = "data/processed/filtered_snps_bpres/lowcov/all.AG.lowcov.merged.filtered.99.20.snps.vcf.gz"
        # Sites with a depth cutoff of 99% percentile and 20% missing data
#        sites = "reports/filtering/depth/lowcov/all.AG.lowcov.{chrom}.filtered.nocall.0.99_0.2.txt",
#        vcf = "data/processed/filtered_snps_bpres/lowcov/all.AG.lowcov.{chrom}.filtered.99.20.snps.vcf.gz"
    output:
#        sites = "data/angsd/lowcov/lowcov.merged.miss20.20k.positions"
        sites = "data/angsd/lowcov/lowcov.{chrom}.miss20.positions"
    shell:
        """
        bcftools query -R {input.sites} -f '%CHROM\t%POS\t%REF\t%ALT{{0}}\n' {input.vcf} > {output.sites}
        """

# (2) Randomly subsample sites sites to 50k sites for the PCA
# shuf -n 50000 data/angsd/lowcov/lowcov.all.miss20.positions >> data/angsd/lowcov/lowcov.all.miss20.50k.positions

# (3) Run the PCA
# -GL 1 calcualtes GLs using GATK method - only calling them so a reference allele can be established for -doMajorMinor
# -doIBS 1 prints a randomly sampled read from each individual at each position
# -doCov 1 prints out the covariance matrix which can be used for a PCA
# -minMaf 0.05 excludes sites with a minor allele freq less the 0.05
# -doMajorMinor 4 specifies the Major allele to be the reference allele specified by -ref

rule PCA_single:
    input:
        ref = config.ref,
        # All Andro
#        bamlist = "data/final_bams/lowcov/all.bamlist"
        # Common garden
        bamlist = "data/final_bams/lowcov/commongarden.bamlist", 
        #  50k random sites, 20% missing data
        sites = "data/angsd/lowcov/lowcov.all.miss20.50k.positions.sorted"
    output:
#        "data/pca/lowcov/all.andro.lowcov.50k.ibs.gz"
        "data/pca/lowcov/cg.andro.lowcov.50k.ibs.gz"
    params:
#        prefix = "data/pca/lowcov/all.andro.lowcov.50k"
        prefix = "data/pca/lowcov/cg.lowcov.50k"
    run:
        shell("angsd \
        -sites {input.sites} \
        -bam {input.bamlist} \
        -doMajorMinor 3 \
        -doCounts 1 \
        -ref {input.ref} \
        -doCov 1 \
        -doIBS 1 \
        -out {params.prefix}")

# (4) Plot in R



###
### Kinship matrix for low coverage common garden ind
###

# (1) extract single reads positions that have all individuals present (no missing data)
# sites file must contain the major and minor (ref and alt) allele
rule cg_ibs:
    input:
        ref = config.ref,
        bamlist = "data/final_bams/lowcov/commongarden.bamlist",
        sites = "data/angsd/lowcov/lowcov.all.miss20.positions"
    output:
        "data/pca/lowcov/cg.andro.lowcov.nomiss.ibs.gz"
    params:
        prefix = "data/pca/lowcov/cg.andro.lowcov.nomiss"
    run:
        shell("angsd \
        -sites {input.sites} \
        -bam {input.bamlist} \
        -minInd 79 \
        -doMajorMinor 3 \
        -doCounts 1 \
        -ref {input.ref} \
        -doIBS 1 \
        -out {params.prefix}")

###
### Kinship matrix for all low coverage
###

###
### Kinship matrix - All Andropogon, low coverage 
###

# (1) Grab a random read at filtered sites with less than 20% missing data for all individuals
rule all_ibs:
    input:
        ref = config.ref,
        # All Andro
        bamlist = "data/final_bams/lowcov/all.bamlist",
        # all sites with less than 20% missing data
        sites = "data/angsd/lowcov/lowcov.all.miss20.positions"
    output:
        "data/pca/lowcov/all.andro.lowcov.all.miss20.ibs.gz"
    params:
        prefix = "data/pca/lowcov/all.andro.lowcov.all.miss20"
    run:
        shell("angsd \
        -sites {input.sites} \
        -bam {input.bamlist} \
        -doMajorMinor 3 \
        -doCounts 1 \
        -ref {input.ref} \
        -doIBS 1 \
        -out {params.prefix}")

# (2) Grab a random read at same sites as step (1) but also require depth > 2
## This needs to be run twice to calculate the diagonal of the kinship matrix
rule kin_diag:
    input:
        ref = config.ref,
        # All Andro
        bamlist = "data/final_bams/lowcov/all.bamlist",
        # all sites with less than 20% missing data
        sites = "data/angsd/lowcov/lowcov.all.miss20.positions"
    output:
        "data/kinship/lowcov/all.andro.lowcov.all.miss20.min6.run{run}.ibs.gz"
    params:
        prefix = "data/kinship/lowcov/all.andro.lowcov.all.miss20.min6.run{run}",
        seed = random.randint(1,100)
    run:
        shell("angsd \
        -sites {input.sites} \
        -bam {input.bamlist} \
        -setMinDepthInd 6 \
        -doMajorMinor 3 \
        -doCounts 1 \
        -ref {input.ref} \
        -doIBS 1 \
        -seed {params.seed} \
        -out {params.prefix}")

# (3) Calculate diagonal (inbreeding) with all sites
rule inbred:
   input:
        "all.andro.lowcov.miss20.min2.100k.run2.ibs.txt"
   output:
        "all.andro.lowcov.miss20.inbreedingcoef.csv"
   run:
        import pandas as pd
        import numpy as np
        run2 = pd.read_table("all.andro.lowcov.miss20.min2.100k.run2.ibs.txt", header = None)
        run1 = pd.read_table("all.andro.lowcov.miss20.min2.100k.run1.ibs.txt", header = None)
        sum = run1 + run2
        sum_mat = sum.iloc[:, 4:]
        af = sum_mat / 2
        af[af < 0] = np.nan
        inbred = af.mean(axis = 0)
        inbred.to_csv("all.andro.lowcov.miss20.inbreedingcoef.csv", index = False, header=False)

###
### STRUCTURE 
###

rule structure:
    input:
        sites = "data/structure/cg.lowcov.50k.structure_input.txt",
        main = "data/structure/mainparams",
        extra = "data/structure/extraparams"
    output:
        "data/structure/cg.lowcov.50k.k{k}.run{run}.75steps.structure_output.txt_f"
    params:
        k = "{k}",
        out = "data/structure/cg.lowcov.50k.k{k}.run{run}.75steps.structure_output.txt"
    run:
        #module load structure-console/2.3.4
        shell("structure -m {input.main} -e {input.extra} -K {params.k} -i {input.sites} -o {params.out}")

###
### Treemix
###

# (2) Convert to treemix format
# treemix_input.R

# (3) Run treemix
rule treemix:
    input:
        "data/treemix/all.andro.lowcov.treemix.txt.gz"
    output:
        "data/treemix/all.andro.lowcov.cov.gz"
    params:
        prefix = "data/treemix/all.andro.lowcov"
    run:
        shell("~/toolsfordayz/treemix-1.13/src/treemix -i {input} -o {params.prefix}")
