### Title: extracting allele depth 
### Author: Alyssa Phillips
### Date: 1/21/2022

library(vcfR, lib.loc = "R_libs")
library("memuse", lib.loc = "R_libs")
library("argparser", lib.loc = "R_libs")

# input is a sorted vcf w/ header for a single ploidy level

# give the arugment a name
ap <- arg_parser("Extracting AD")

# add mandatory positional arguments (filename)
ap <- add_argument(ap, "vcf", help = "Input vcf - sorted vcf w/ header containg genotypes of single ploidy")

# additional arguments
ap <- add_argument(ap, "--chr", help = "scaffold number")
ap <- add_argument(ap, "--ploidy", help = "ploidy")

# parse arguments
argv <- parse_args(ap)

# Input
vcf_file <- as.character(argv$vcf)
chr <- as.character(argv$chr)
p <- as.character(argv$p)

vcf_file 

#error <- as.numeric(snakemake@input[[2]])
error <- 0.01

# load vcf
vcf <- read.vcfR(vcf_file)

print("loaded vcf")

# grab genotype depth
dp <- extract.gt(vcf, element = 'DP', as.numeric = T) # extract as numeric

print("extracted DP")

# grab allele depth
ad <- extract.gt(vcf, element = 'AD', as.numeric = F)

print("grab AD")

# extract alt (record = 2) read counts
alt <- masplit(ad,  record = 2, sort = 0)

print("split AD")

# transform to be loci = columns & genotypes = rows
t_dp <- t(dp)
t_alt <- t(alt)

# replace zeros (missing data) with -9 in the total read count matrix
t_dp[t_dp == 0] <- -9

# replace sites with missing data as -9 while leaving true zeros alone
t_alt[t_dp == -9] <- -9

print("transformed and replace NA with -9")
# write per locus error rate vector
E <- matrix(data = error, nrow = dim(t_alt)[2], ncol = 1)

# export sample names and SNP positions (so we know their order later)
geno <- rownames(t_alt)
snps <- colnames(t_alt)

# Export data
getwd()
paste0("data/ebg/lowcov/total_reads.",chr, ".", p, ".txt")
write.table(t_dp, file = paste0("./data/ebg/lowcov/total_reads.",chr, ".", p, ".txt"), sep = "\t", row.names = F, col.names = F)
print("wrote total reads mat")
write.table(t_alt, file = paste0("data/ebg/lowcov/alt_reads.",chr, ".", p, ".txt"), sep = "\t", row.names = F, col.names = F)
print("write alt mat")
write.table(E, file = paste0("data/ebg/lowcov/error_rate.",chr, ".", p, ".txt"), sep = "\t", row.names = F, col.names = F)
write.table(geno, file = paste0("data/ebg/lowcov/samples.",chr, ".", p, ".txt"), sep = "\t", row.names = F, col.names = F)
write.table(snps, file = paste0("data/ebg/lowcov/snp_positions.",chr, ".", p, ".txt"), sep = "\t", row.names = F, col.names = F)
