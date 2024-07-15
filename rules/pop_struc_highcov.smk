# Analyses on 6x high coverage individuals

###
### PCA
###

# Run for verification that single-read PCA is correct

# (1) Estimate GL in beagle format
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

# (2) Merge beagle files together
#        shell("cat <(zcat data/angsd/pca/*beagle.gz | head -n1) <(zcat data/angsd/pca/*beagle.gz | \
#        grep -v marker) > data/angsd/pca/highcov.merged.8dp70.beagle
#        shell("gzip data/angsd/pca/highcov.merged.8dp70.beagle")

# (3) Grab 30k random snps
#        zcat highcov.merged.8dp70.beagle.gz | head -n1 > highcov.merged.8dp70.30k.beagle
#        shuf -n 30000 <(zcat highcov.merged.8dp70.beagle.gz) >> highcov.merged.8dp70.30k.beagle
#        gzip highcov.merged.8dp70.30k.beagle

# (4) Run PCA
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
### Estimate the inbreeding coefficient (Fis)
###

# (1) Call genotype likelihoods with ANGSD

# -doGlf 3 required for ngsF - binary 3 times likelihood (.glf.gz)
# -GL 1 - SAMtools GL model
# -doMaf 1 - minor allele frequency is known
# -doMajorMinor 4 - known major allele is from the reference
# Don't specify a minimum number of individuals as it is not necessary

rule call_gls:
    input:
#        bams = expand("data/final_bams/6x_subsample/{low_geno}_{low_per}.subsample.bam", zip, low_geno = LOW_GEN>
        ref = config.ref,
        bamlist = "data/final_bams/highcov/highcov.6x.bamlist"
    output:
        "data/angsd/highcov_6x/highcov_6x_andro.{chrom}.glf.gz"
    params:
        prefix = "data/angsd/highcov_6x/highcov_6x_andro.{chrom}",
        chrom = "{chrom}"
    shell:
        """
        angsd \
        -GL 1 -P 15 \
        -doGlf 3 \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
        -minMapQ 30 -minQ 30 \
        -doCounts 1 -setMinDepthInd 8 -setMaxDepthInd 70 \
        -b {input.bamlist} \
        -r {params.chrom} \
        -doMaf 1 \
        -doMajorMinor 4 \
        -ref {input.ref} \
        -SNP_pval 1e-6 \
        -out {params.prefix}
        """

# (2) Estimate inbreeding coefficient (F) with ngsF

## Use approximated method first and then main EM algorithm if needed
## Should be run multiple times to assess convergence

# --init_values r - start with random intial values (recommended for low coverage)
# --min_epsilon - Maximum RMSD between iterations to assume convergence (ideally want a value of zero)
rule ngsF:
    input:
        gl = "data/angsd/highcov_6x/highcov_6x_andro.{chrom}.glf.gz",
        maf = "data/angsd/highcov_6x/highcov_6x_andro.{chrom}.mafs.gz"
    output:
        "data/angsd/highcov_6x/highcov_6x_andro.{chrom}.run{run}.approx_indF"
    shell:
        """
        NSITES=`zcat {input.maf} | tail -n+2 | wc -l`
        echo $NSITES

        SEED=$RANDOM
        echo $SEED

        zcat {input.gl} | \
        ngsF \
        --n_ind 47 \
        --n_sites $NSITES \
        --init_values r \
        --glf - \
        --out {output} \
        --approx_EM \
        --seed $SEED \
        --min_epsilon 1e-6 \
        --n_threads 52
        """

# (3) Run deep search with estimated starting values
rule ngsF_deep:
    input:
        gl = "data/angsd/highcov_6x/highcov_6x_andro.{chrom_avail}.glf.gz",
        init = "data/angsd/highcov_6x/highcov_6x_andro.{chrom_avail}.run{run}.approx_indF.pars",
        maf = "data/angsd/highcov_6x/highcov_6x_andro.{chrom_avail}.mafs.gz"
    output:
        "data/angsd/highcov_6x/highcov_6x_andro.{chrom_avail}.run{run}.indF"
    shell:
        """
        NSITES=`zcat {input.maf} | tail -n+2 | wc -l`
        echo $NSITES

        SEED=$RANDOM
        echo $SEED

        zcat {input.gl} | \
        ngsF \
        --n_ind 47 \
        --n_sites $NSITES \
        --init_values r \
        --glf - \
        --out {output} \
        --init_values {input.init} \
        --seed $SEED \
        --min_epsilon 1e-6 \
        --n_threads 52 \
        """

# (4) Plot in R - bigblue_Fis.R
