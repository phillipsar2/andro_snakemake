import config
import pandas as pd

# Different SAMPLE names
#SAMPLE = glob_wildcards("/group/jrigrp10/andropogon_shortreads/{sample}.merge.R1.fastq.gz").sample
#SAMPLE = glob_wildcards("/group/jrigrp10/andropogon_shortreads/{sample}_1.fq.gz").sample
SAMPLE = glob_wildcards("/group/jrigrp10/andropogon_shortreads/ucd_seq/{sample}_R1_001.fastq.gz").sample
#print(SAMPLE)

# BAM names encompass all samples regardless of their fastq pattern
BAM = glob_wildcards("data/interm/mark_dups/{bam}.dedup.bam").bam
#BAM = glob_wildcards("data/final_bams/lowcov/{bam}.dedup.bam").bam
#BAM = glob_wildcards("data/final_bams/highcov/{bam}.dedup.bam").bam
#print(BAM)

# MERGE contains a list of the bams that belong to each genotype (GENO) so they can be merged
file = pd.read_csv("bams_to_merge.tsv", sep = "\t", header = 0)
MERGE_A = list(file.Merge_A)
MERGE_B = list(file.Merge_B)
GENO = list(file.Genotype)

# HIGH contains the bams to be subsampled to low coverage
HIGH = glob_wildcards("data/interm/mark_dups/IN{high}.dedup.bam").high

# LOW_SUB contains the lowcov processed bams that need to be subsampled to even lower coverage


# filep contains the bams nQuire was run on to determine ploidy
filep = pd.read_csv("bams_nquire.csv", sep = "\t", header = None)
UNKNOWN = list(filep[0])

# List of genotypes and ploidy levels
#fileg = pd.read_csv(config.highcov, header = 0)
#PLOIDY = list(fileg.Ploidy)
#VCF = list(fileg.Sequencefile)

# List of all chromsomes
#chr = pd.read_csv(config.contig_list, header = None)
#CHROM = list(chr[0])
## low cov scaffolds w sites not filtered out (some lost entirely)
CHROM = ["Chr01A","Chr01B","Chr01C","Chr02A","Chr02B","Chr02C","Chr03A","Chr03B","Chr03C","Chr04A","Chr04B","Chr04C","Chr05A","Chr05B","Chr05C","Chr06A","Chr06B","Chr06B","Chr06C","Chr07A","Chr07B","Chr07C","Chr08A","Chr08C","Chr08B","Chr09A","Chr09B","Chr09C","Chr10A","Chr10B","Chr10C","scaffold_144","scaffold_163","scaffold_32","scaffold_490","scaffold_542","scaffold_965"]
## high cov scaffolds w sites not filtered out (some lost entirely)
#CHROM =  ["Chr01A", "Chr01B", "Chr01C", "Chr02A", "Chr02C", "Chr02B", "Chr03A", "Chr03B", "Chr03C", "Chr04A", "Chr04B", "Chr04C", "Chr05A", "Chr05B", "Chr05C", "Chr06A", "Chr06B", "Chr06C", "Chr07A", "Chr07B", "Chr07C", "Chr08A", "Chr08B", "Chr08C", "Chr09A", "Chr09B", "Chr09C", "Chr10A", "Chr10B", "Chr10C", "scaffold_144", "scaffold_163", "scaffold_263", "scaffold_32", "scaffold_490", "scaffold_542", "scaffold_590", "scaffold_64", "scaffold_965"]
#print(CHROM)

# Select working with high or low coverge samples for SNP calling 
# options = ["highcov", "lowcov"]
COV = ["lowcov"]

# Set SNP filtering parameters
p = ["99"]
miss = ["20"]

# Set what ploidy group working with for GL calling
# options = ["9x", "6x"]
CYT = ["6x"]

# number of admixture groups (k) - for 14 groups do range(15))
K = list(range(2,15))

# which chain are we currently running for entropy?
CHAIN = ["1"]

# Rule all describes the final output of the pipeline
rule all:
    input:
        ## Trimming raw reads
#         expand("reports/fastp/{sample}.json", sample = SAMPLE)
        ## Aligning reads
#        expand("data/interm/mark_dups/{sample}.dedup.bam", sample = SAMPLE),
#        expand("reports/bamqc/{bam}_stats/qualimapReport.html", bam = BAM),
        ## Merging low-coverage bams
#        expand("data/interm/mark_dups/{geno}.{merge_A}.{merge_B}.merged.dedup.bam", zip, merge_A = MERGE_A, merge_B = MERGE_B, geno = GENO),
#        expand("data/interm/mark_dups/{geno}.{merge_A}.{merge_B}.merged.rg.dedup.bam", zip, merge_A = MERGE_A, merge_B = MERGE_B, geno = GENO),
#        expand("reports/bamqc/merged/{geno}.{merge_A}.{merge_B}_stats/genome_results.txt",zip, merge_A = MERGE_A, merge_B = MERGE_B, geno = GENO)
        ## Subsample high-coverage bams
#        expand("data/final_bams/lowcov/IN{high}.subsample.dedup.bam", high = HIGH),
#        expand("reports/bamqc/final_bams/lowcov/{bam}_stats/genome_results.txt", bam = BAM),
        ## nQuire
#        hist = expand("data/nQuire/{unknown}_denoised.hist", unknown = UNKNOWN),
#        r = "data/nQuire/nquire_results_denoised.txt",
        ## SNP Calling
#        index = expand("data/final_bams/lowcov/{bam}.dedup.bam.bai", bam = BAM),
#        mpileup = expand("data/vcf/highcov/all.AG.highcov.{chr}.raw.vcf.gz", chr = CHROM)
#        mpileup = expand("data/vcf/lowcov/ucd.lowcov.{chr}.raw.vcf.gz", chr = CHROM)
        ## Filtering
#        index = expand("data/vcf/{cov}/all.AG.{cov}.{chr}.raw.vcf.gz.tbi", chr = CHROM, cov = COV),
#        snp = expand("data/raw/vcf_bpres/{cov}/all.AG.{cov}.{chr}.raw.snps.vcf.gz", chr = CHROM, cov = COV),
#        diag = expand("reports/filtering/{cov}/all.AG.{cov}.{chr}.table", chr = CHROM, cov = COV)
#        hf = expand("data/processed/filtered_snps_bpres/{cov}/all.AG.{cov}.{chr}.filtered.nocall.vcf", chr = CHROM, cov = COV),
#        dp_diag = expand("reports/filtering/depth/{cov}/all.AG.{cov}.{chr}.filtered.nocall.table", chr = CHROM, cov = COV),
#        dp_miss = expand("reports/filtering/depth/{cov}/all.AG.{cov}.{chr}.filtered.nocall.0.99_0.2.txt", p = p, miss = miss, chr = CHROM, cov = COV)
#        dp_miss = expand("reports/filtering/depth/{cov}/all.AG.{cov}.{chr}.filtered.nocall.0.99_0.2_8.txt", chr = CHROM, cov = COV),
         dp_miss = expand("reports/filtering/depth/{cov}/all.AG.{cov}.{chr}.filtered.nocall.0.99_0_1.txt", chr = CHROM, cov = COV),
#        grab_snps = expand("data/processed/filtered_snps_bpres/{cov}/all.AG.{cov}.{chr}.filtered.{p}.{miss}.snps.vcf.gz", cov = COV, p = p, miss = miss, chr = CHROM),
#        grab_gl = expand("data/processed/filtered_snps_bpres/{cov}/all.AG.{cov}.{chr}.PL.txt", cov = COV, chr = CHROM)
        ## EBG
#        split_ploidy = expand("data/processed/filtered_snps_bpres/lowcov/AG.lowcov.{chr}.{ploidy}.snps.vcf", chr = CHROM, ploidy = CYT),
#        ad_mat = expand("data/ebg/lowcov/total_reads.{chr}.{ploidy}.txt", chr = CHROM, ploidy = CYT),
#        ebg = expand("data/ebg/lowcov/{chr}.{ploidy}-PL.txt", chr = CHROM, ploidy = CYT),
#        gl_mat = expand("data/ebg/lowcov/{chr}-GL.txt", chr = CHROM),
#        rand_snps = "data/ebg/lowcov/10k_lowcov-GL.txt"
#        filt_miss = expand("data/ebg/lowcov/genoliks/{chr}.miss20-GL.txt", chr = CHROM)
        ## ENTROPY
#         ent_in = "data/entropy/lowcov/10k_lowcov.miss5.mpgl",
#        ent_in_tab = "data/entropy/highov/10k_highcov.mpgl"
#        entropy = expand("data/entropy/{cov}/10k_{cov}.{k}.{chain}.hdf5", cov = COV, k = K, chain = CHAIN)
#        ent_out = expand("data/entropy/{cov}/qest.k{k}.c{chain}.txt", cov = COV, k = K, chain = CHAIN)
        ## High cov PCA
#        pca =  expand("data/angsd/pca/highcov.{chrom}.8dp70.beagle.gz", chrom = CHROM)
#        pcaangsd = "data/angsd/pca/highcov.merged.8dp70.cov"
        ## Low cov PCA
#        get_sites = expand("data/angsd/lowcov/lowcov.{chrom}.miss20.positions", chrom = CHROM)
#        all_single = "data/angsd/lowcov/all.andro.lowcov.all.positions.ibs.gz",
#        pca_single = "data/pca/lowcov/all.andro.lowcov.50k.ibs.gz"
#        pca_single = "data/pca/lowcov/cg.andro.lowcov.50k.ibs.gz"
        ## Kinship matrix
         cg_ibs = "data/pca/lowcov/cg.andro.lowcov.nomiss.ibs.gz"

## Rules
#include: "rules/mapping.smk"
#include: "rules/process_bam.smk"
#include: "rules/determine_ploidy.smk"
#include: "rules/calling.smk"
#include: "rules/filtering.smk"
#include: "rules/gl_calling.smk"
include: "rules/pop_struc.smk"
