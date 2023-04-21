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
        ibs = "data/kinship/lowcov/all.andro.lowcov.all.miss20.min2.run{run}.ibs.gz",
        bamrand = "data/kinship/lowcov/bamlist.rand.run{run}.txt"
    params:
        prefix = "data/kinship/lowcov/all.andro.lowcov.all.miss20.min2.run{run}",
        seed = random.randint(1,100)
    run:
        shell("shuf {input.bamlist} > {output.bamrand}")
        shell("angsd \
        -sites {input.sites} \
        -bam {output.bamrand} \
        -setMinDepthInd 2 \
        -doMajorMinor 3 \
        -doCounts 1 \
        -ref {input.ref} \
        -doIBS 1 \
        -seed {params.seed} \
        -out {params.prefix}")

# (3) Subset random 100k SNPs
rule kin_sub:
    input:
        run1 = "data/kinship/lowcov/all.andro.lowcov.all.miss20.min2.run1.ibs.gz",
        run2 = "data/kinship/lowcov/all.andro.lowcov.all.miss20.min2.run2.ibs.gz"
    output:
        run1 = "data/kinship/lowcov/all.andro.lowcov.miss20.min2.100k.run1.ibs.txt",
        run2 = "data/kinship/lowcov/all.andro.lowcov.miss20.min2.100k.run2.ibs.txt"
    shell:
        """
        get_seeded_random(){{ \\
        seed="$1" \\
        openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \\
        </dev/zero 2>/dev/null}}
        shuf -n 100000 --random-source=<(get_seeded_random 24) <(zcat {input.run1) > {output.run1}
        shuf -n 100000 --random-source=<(get_seeded_random 24) <(zcat {input.run2) > {output.run2} 
        """


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
#        sites = "data/structure/cg.lowcov.50k.structure_input.txt",
        sites = "data/structure/all.andro.lowcov.100k.structure_input.txt",
        main = "data/structure/mainparams",
        extra = "data/structure/extraparams"
    output:
#        "data/structure/cg.lowcov.50k.k{k}.run{run}.75steps.structure_output.txt_f"
        "data/structure/al.andro.lowcov.100k.k{k}.run{run}.75steps.structure_output.txt_f"
    params:
        k = "{k}",
#        out = "data/structure/cg.lowcov.50k.k{k}.run{run}.75steps.structure_output.txt"
        out = "data/structure/all.andro.lowcov.100k.k{k}.run{run}.75steps.structure_output.txt"
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
