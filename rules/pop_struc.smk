import random

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

# (4) Plot in R - scripts local


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

# (2) Estimate kinship matrix in custom script in R

###
### STRUCTURE 
###

rule structure:
    input:
        ## Common garden
#        sites = "data/structure/cg.lowcov.50k.structure_input.txt",
        ## All andropogon
#        sites = "data/structure/all.andro.lowcov.100k.structure_input.txt",
        ## All andropogon without inbreds
        sites = "data/structure/noinbreds.andro.lowcov.100k.structure_input.txt",
        main = "data/structure/mainparams",
        extra = "data/structure/extraparams"
    output:
#        "data/structure/cg.lowcov.50k.k{k}.run{run}.75steps.structure_output.txt_f"
#        "data/structure/all.andro.lowcov.100k.k{k}.run{run}.75steps.structure_output.txt_f"
        "data/structure/noinbreds.andro.lowcov.100k.k{k}.run{run}.75steps.structure_output.txt_f"
    params:
        k = "{k}",
#        out = "data/structure/cg.lowcov.50k.k{k}.run{run}.75steps.structure_output.txt"
#        out = "data/structure/all.andro.lowcov.100k.k{k}.run{run}.75steps.structure_output.txt"
        out = "data/structure/noinbreds.andro.lowcov.100k.k{k}.run{run}.75steps.structure_output.txt"
    run:
        #module load structure-console/2.3.4
        shell("structure -m {input.main} -e {input.extra} -K {params.k} -i {input.sites} -o {params.out}")

###
### Inbreeding in 6x
###

# (1) call genotype likelihoods with ANGSD

# -doGlf 3 required for ngsF - binary 3 times likelihood (.glf.gz)
# -GL 1 - SAMtools GL model
# -doMaf 1 - minor allele frequency is known
# -doMajorMinor 4 - known major allele is from the reference
# Don't specify a minimum number of individuals as it is not necessary

rule call_gls:
    input:
#        bams = expand("data/final_bams/6x_subsample/{low_geno}_{low_per}.subsample.bam", zip, low_geno = LOW_GENO, low_per = LOW_PER),
        ref = config.ref,
        bamlist = "data/final_bams/6x_subsample/6x_subsampled.bamlist"
    output:
        "data/angsd/lowcov_6x/lowcov_6x_andro.{chrom}.glf.gz"
    params:
        prefix = "data/angsd/lowcov_6x/lowcov_6x_andro.{chrom}",
        chrom = "{chrom}"
    shell:
        """
        angsd \
        -GL 1 -P 15 \
        -doGlf 3 \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
        -minMapQ 30 -minQ 30 \
        -doCounts 1 -setMinDepthInd 1 -setMaxDepthInd 6 \
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
        gl = "data/angsd/lowcov_6x/lowcov_6x_andro.all.glf"
    output:
        "data/angsd/lowcov_6x/lowcov_6x_andro.run1.approx_indF"
    shell:
        """
        ngsF \
        --n_ind 107 \
        --n_sites 45357448 \
        --init_values r \
        --glf {input.gl} --out {output} \
        --approx_EM \
        --seed 12345 \
        --min_epsilon 1e-6
        """

# (3) Run deep search with estimated starting values
rule ngsF_deep:
    input:
        gl = "data/angsd/lowcov_6x/lowcov_6x_andro.all.glf",
        init = "data/angsd/lowcov_6x/lowcov_6x_andro.run1.approx_indF.pars"
    output:
        "data/angsd/lowcov_6x/lowcov_6x_andro.run1.indF"
    shell:
        """
        ngsF \
        --n_ind 107 \
        --n_sites 45357448 \
        --init_values r \
        --glf {input.gl} --out {output} \
        --init_values {input.init} \
        --seed 1235 \
        --min_epsilon 1e-6 \
        --n_threads 30 \
        """
    
    
###
### Diversity stats
###

# (1) Estimate the global SFS 
# -doSaf 1 estimate the SFS based on ind genotype likelihoods
# -minInd 3 100% of the genotypes (n = 3) need data at the site
rule angsd_saf:
    input:
        ref = config.ref,
        bamlist = "data/final_bams/6x_subsample/bamlists/{pops}.bamlist"
    output:
        temp("data/angsd/saf/{pops}.{chrom}.50per.saf.gz")
    params:
        prefix = "data/angsd/saf/{pops}.{chrom}.50per",
        chrom = "{chrom}"
    shell:
        """
        angsd \
        -GL 1 -P 15 \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
        -minMapQ 30 -minQ 30 \
        -doCounts 1 -setMinDepthInd 1 -setMaxDepthInd 6 \
        -b {input.bamlist} \
        -r {params.chrom} \
        -minInd 3 \
        -ref {input.ref} -anc {input.ref} \
        -doSaf 1 \
        -out {params.prefix}
        """

## (2) Generate global estimate of SFS (site frequency spectrum) then per-site thetas
# -fold 1 specifies the folded spectrum as I don't have an ancestral state
# real SFS saf2theta calculates the thetas for each site
rule pop_sfs:
    input:
        ref = config.ref,
        saf = "data/angsd/saf/{pops}.{chrom}.50per.saf.gz",
    output:
        sfs = temp("data/angsd/saf/{pops}.{chrom}.sfs")
    params:
        prefix = "data/angsd/saf/{pops}.{chrom}"
    run:
        shell("realSFS {params.prefix}.saf.idx  -P 10 -fold 1 > {params.prefix}.sfs")
        shell("realSFS saf2theta {params.prefix}.saf.idx -sfs {params.prefix}.sfs -outname {params.prefix}")

# (3) calculate thetas (and neutrality tests) in 10k sliding windows
### when using a folded SFS, only thetaW (tW), thetaD (tP), and tajimasD will be meaningful in the output of realSFS
rule pop_pi:
    input:
        sfs = "data/angsd/saf/{pops}.{chrom}.sfs"
    output:
        stats = "data/angsd/theta/{pops}.{chrom}.{window}.thetas.idx.pestPG"
    params:
        prefix = "data/angsd/theta/{pops}.{chrom}.{window}",
        win = "{window}"
    run:
        shell("thetaStat do_stat {params.prefix}.thetas.idx \
        -win {params.win} \
        -step {params.win}")
       # merge the files
#        shell("cat <(cat *.thetas.idx.pestPG | head -n1) <(cat *.thetas.idx.pestPG | grep -v nSites) > all.boulder.thetas.idx.pestPG")
