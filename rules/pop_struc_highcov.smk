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


