### Title: extracting allele depth 
### Author: Alyssa Phillips
### Date: 1/21/2022

library(vcfR)

# input is a sorted vcf w/ header for a single ploidy level

# Snakemake syntax inputs
vcf_file <- as.character(snakemake@input[[1]])
error <- as.numeric(snakemake@input[[2]])

# load vcf
vcf <- read.vcfR(vcf_file)

# grab genotype depth
dp <- extract.gt(vcf, element = 'DP', as.numeric = T) # extract as numeric

# grab allele depth
ad <- extract.gt(vcf, element = 'AD', as.numeric = F)

# extract alt (record = 2) read counts
alt <- masplit(ad,  record = 2, sort = 0)

# transform to be loci = columns & genotypes = rows
t_dp <- t(dp)
t_alt <- t(alt)

# replace zeros (missing data) with -9 in the total read count matrix
t_dp[t_dp == 0] <- -9

# replace sites with missing data as -9 while leaving true zeros alone
t_alt[t_dp == -9] <- -9

# write per locus error rate vector
E <- matrix(data = error, nrow = dim(t_alt)[2], ncol = 1)

# export sample names and SNP positions (so we know their order later)
geno <- rownames(t_alt)
snps <- colnames(t_alt)

# export in snakemake format
write.table(t_dp, file = snakemake@output[[1]], sep = "\t", row.names = F, col.names = F)
write.table(t_alt, file = snakemake@output[[2]], sep = "\t", row.names = F, col.names = F)
write.table(E, file = snakemake@output[[3]], sep = "\t", row.names = F, col.names = F)
write.table(geno, file = snakemake@output[[4]], sep = "\t", row.names = F, col.names = F)
write.table(snps, file = snakemake@output[[5]], sep = "\t", row.names = F, col.names = F)
