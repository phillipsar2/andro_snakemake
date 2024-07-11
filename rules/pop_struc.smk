import random

###
### PCA on low coverage dataset - single read sampling approach
###

# (1) Create the sites files
rule create_sites:
    input:
        # Sites with a depth cutoff of 99% percentile and 20% missing data
#        sites = "reports/filtering/depth/lowcov/all.AG.lowcov.{chrom}.filtered.nocall.0.99_0.2.txt",
#        vcf = "data/processed/filtered_snps_bpres/lowcov/all.AG.lowcov.{chrom}.filtered.99.20.snps.vcf.gz"
        # No clones, 20% missing data and 99% percentile
        sites = "reports/filtering/depth/{cov}/all.AG.noclones.{cov}.{chr}.filtered.nocall.0.99_0.2.txt",
        vcf = "data/processed/filtered_snps_bpres/{cov}/all.AG.noclones.{cov}.{chr}.filtered.nocall.vcf"
    output:
#        sites = "data/angsd/lowcov/lowcov.{chrom}.miss20.positions"
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
         bamlist = "data/final_bams/lowcov/all.noclones.bamlist",
        # Common garden
#        bamlist = "data/final_bams/lowcov/commongarden.bamlist", 
        #  50k random sites, 20% missing data
#        sites = "data/angsd/lowcov/lowcov.all.miss20.50k.positions.sorted"
        # 100k random sites, 20% missing data cutoff - 3/10/2024
         sites = "data/angsd/lowcov/lowcov.all.miss20.100k.sorted.positions"
        # All sites with 20% missing data cutoff - 3/7/2024
#         sites = "data/angsd/lowcov/lowcov.all.miss20.positions"
    output:
#        "data/pca/lowcov/all.andro.lowcov.50k.ibs.gz"
#        "data/pca/lowcov/cg.andro.lowcov.50k.ibs.gz"
         "data/pca/lowcov/all.andro.lowcov.miss20.100k.covMat"
    params:
#        prefix = "data/pca/lowcov/all.andro.lowcov.50k"
#        prefix = "data/pca/lowcov/cg.lowcov.50k"
         prefix = "data/pca/lowcov/all.andro.lowcov.miss20.100k"
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
#        sites = "data/structure/noinbreds.andro.lowcov.100k.structure_input.txt",
        ## All andropogon without clones
        sites = "data/structure/all.andro.noclones.lowcov.100k.2024-03-08.structure_input.txt",
        main = "data/structure/mainparams",
        extra = "data/structure/extraparams"
    output:
#        "data/structure/cg.lowcov.50k.k{k}.run{run}.75steps.structure_output.txt_f"
#        "data/structure/all.andro.lowcov.100k.k{k}.run{run}.75steps.structure_output.txt_f"
#        "data/structure/noinbreds.andro.lowcov.100k.k{k}.run{run}.75steps.structure_output.txt_f"
        "data/structure/all.andro.noclones.lowcov.100k.k{k}.run{run}.75steps.structure_output.txt_f"
    params:
        k = "{k}",
#        out = "data/structure/cg.lowcov.50k.k{k}.run{run}.75steps.structure_output.txt"
#        out = "data/structure/all.andro.lowcov.100k.k{k}.run{run}.75steps.structure_output.txt"
#        out = "data/structure/noinbreds.andro.lowcov.100k.k{k}.run{run}.75steps.structure_output.txt"
        out = "data/structure/all.andro.noclones.lowcov.100k.k{k}.run{run}.75steps.structure_output.txt"
    run:
        #module load structure-console/2.3.4
        shell("structure -m {input.main} -e {input.extra} -K {params.k} -i {input.sites} -o {params.out}")

###
### Thetas - 6x only
###

# (1) Estimate the SFS for each population
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
        sfs = "data/angsd/saf/{pops}.{chrom}.sfs"
    params:
        prefix = "{pops}.{chrom}"
    shell:
        """
        realSFS data/angsd/saf/{params.prefix}.50per.saf.idx  -P 10 -fold 1 > data/angsd/saf/{params.prefix}.sfs
#        realSFS saf2theta data/angsd/saf/{params.prefix}.saf.idx -sfs {params.prefix}.sfs -outname data/angsd/saf/{params.prefix}
        """

# (3) calculate thetas (and neutrality tests) in 10k sliding windows
### when using a folded SFS, only thetaW (tW), thetaD (tP), and tajimasD will be meaningful in the output of realSFS
rule pop_pi:
    input:
        sfs = "data/angsd/saf/{pops}.{chrom}.sfs",
        saf = "data/angsd/saf/{pops}.{chrom}.50per.saf.idx"
    output:
        stats = "data/angsd/saf/{pops}.{chrom}.{window}.thetas.idx.pestPG"
    params:
        prefix = "data/angsd/saf/{pops}.{chrom}.{window}",
        win = "{window}"
    shell:
        """
        realSFS saf2theta {input.saf} -sfs {input.sfs} -outname {params.prefix}
        if test -s {params.prefix}.thetas.idx; \
        then \
            thetaStat do_stat {params.prefix}.thetas.idx -win {params.win} -step {params.win}; \
        else \
            touch {output}; \
        fi
#        thetaStat do_stat {params.prefix}.thetas.idx \
#        -win {params.win} \
#        -step {params.win}
#        if test -f {params.prefix}.thetas.idx; then echo "idx exists"; else touch {output}; fi
        """
       # merge the files
#        shell("cat <(cat *.thetas.idx.pestPG | head -n1) <(cat *.thetas.idx.pestPG | grep -v nSites) > all.boulder.thetas.idx.pestPG")


## (4) Merge the theta estimates
rule merge_pi:
    input:
        expand("data/angsd/saf/{pops}.{chrom}.{window}.thetas.idx.pestPG", pops = POPS, chrom = CHROM, window = WINDOW)
    output:
        "data/angsd/saf/{pops}.all.{window}.thetas.idx.pestPG.gz"
    params:
        "data/angsd/saf/{pops}.all.{window}.thetas.idx.pestPG"
    shell:
        """
        cat <(cat data/angsd/saf/{wildcards.pops}.*.{wildcards.window}.thetas.idx.pestPG| head -n1 | cut -f 1,2,4,5,9,14 ) \
        <(cat data/angsd/saf/{wildcards.pops}.*.{wildcards.window}.thetas.idx.pestPG | grep -v nSites | cut -f 1,2,4,5,9,14 | awk '$6 != 0') | \
        gzip > {output}
        """

###
### Heterozygosity - 6x only
###

## (1) Estimate the folded SFS for all hexaploids range-wide
### Heterozygosity is the second bin of the SFS
# -minInd 106 is 20% missing data per site
# -fold 1 estimate folded SFS
rule all_saf:
    input:
        ref = config.ref,
        bamlist = "data/final_bams/6x_subsample/bamlists/all.6x.bamlist"
    output:
        "data/angsd/saf/all.6x.{chrom}.saf.gz"
    params:
        prefix = "data/angsd/saf/{pops}.{chrom}.",
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
        -minInd 106
        -ref {input.ref} -anc {input.ref} \
        -doSaf 1 \
        -fold 1 \
        -out {params.prefix}
        """

rule all_sfs:
     input:
         "data/angsd/saf/all.6x.{chrom}.saf.idx"
     output:
         "data/angsd/saf/all.6x.{chrom}.est.ml"
     shell:
         """
         realSFS {input} > {output}
         """
