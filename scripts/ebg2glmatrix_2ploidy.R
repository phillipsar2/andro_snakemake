### Title: Convert EBG output format to a GL matrix - two ploidy levels
### Author: Alyssa Phillips
### Date: 3/18/22

library(stringr)
library(MASS)
library(dplyr)
library("argparser", lib.loc = "R_libs")

# > commandArgs ----
# give the arugment a name
ap <- arg_parser("convert EBG PL file to GL mat - two ploidy levels")

# add mandatory positional arguments (filename)
ap <- add_argument(ap, "pl", help = "PL file output from EBG for first ploidy level")
ap <- add_argument(ap, "pl2", help = "PL file output from EBG for second ploidy level")

# additional arguments
ap <- add_argument(ap, "--tot", help = "total read depth file for first ploidy level that is the input for EBG")
ap <- add_argument(ap, "--ploidy", help = "ploidy of first sample file; single numeric value")
ap <- add_argument(ap, "--geno", help = "list of samples names for first ploidy level  output from vcf2ADmatrix.R script")
ap <- add_argument(ap, "--snps", help = "list of SNP positions for either ploidy level output from vcf2ADmatrix.R script - should be the same for both files")
ap <- add_argument(ap, "--tot2", help = "total read depth file for second ploidy level that is the input for EBG")
ap <- add_argument(ap, "--ploidy2", help = "ploidy of second sample file; single numeric value")
ap <- add_argument(ap, "--geno2", help = "list of samples names for the second ploidy level output from vcf2ADmatrix.R script")

# parse arguments
argv <- parse_args(ap)


# > Input data ----
# adjust to be arguments from the command (need to learn how to do this)
# pl_file <- as.character("~/Andropogon/EBG/Chr07B.6x.10k-PL.txt")
# reads_file <- as.character("~/Andropogon/EBG/total_reads.Chr07B.6x.10k.txt")
# p = 2
# 
# pl_file2 <- as.character("~/Andropogon/EBG/Chr07B.9x.10k-PL.txt")
# reads_file2 <- as.character("~/Andropogon/EBG/total_reads.Chr07B.9x.10k.txt")
# p2 = 3

# # metadata
# genos <- read.csv("~/Andropogon/EBG/samples.Chr07B.6x.10k.txt", sep = "\t", header = F)
# pos <- read.csv("~/Andropogon/EBG/snp_positions.Chr07B.6x.10k.txt", sep = "\t", header = F)
# 
# genos2 <- read.csv("~/Andropogon/EBG/samples.Chr07B.9x.10k.txt", sep = "\t", header = F)
# #pos2 <- read.csv("~/Andropogon/EBG/snp_positions.Chr07B.9x.10k.txt", sep = "\t", header = F)

pl_file <- as.character(argv$pl)
reads_file <- as.character(argv$tot)
p = as.numeric(argv$ploidy)


pl_file2 <- as.character(argv$pl2)
reads_file2 <- as.character(argv$tot2)
p2 = as.numeric(argv$ploidy2)

fname <- strsplit(pl_file, ".6x-PL.txt")[[1]]
fname2 <- strsplit(pl_file2, ".9x-PL.txt")[[1]]

# metadata
gfile <- as.character(argv$geno)
genos <- read.csv(gfile, sep = "\t", header = F)
pfile <- as.character(argv$s)
pos <- read.csv(pfile, sep = "\t", header = F)
head(pos)

gfile2 <- as.character(argv$geno2)
genos2 <- read.csv(gfile2, sep = "\t", header = F)
head(genos2)

print("read in metadata")

# >  Read in first PL file line by line and split columns by tab separation ----
pl <- readLines(pl_file) %>% str_split(pattern = "\t", simplify = TRUE) 
pl <- apply(pl, 2, function(x) gsub("^$|^ $", NA, x)) # replace blanks with NA

print("read in pl")

# load EBG input -- the total read counts matrix
read_counts_file <- read.csv(reads_file, sep = "\t", header = F )
read_counts <- t(read_counts_file)

print("read in readcounts")

# > Edit first file ----
# If values in row != -9,  assign PL to those spots in the new array [R, C] = [loci, individuals] 
pl_mat <- matrix(data = NA,
                 nrow = dim(read_counts)[1],
                 ncol = dim(read_counts)[2] ) # make an empty matrix

l = dim(read_counts)[1]
for (i in 1:l){
  pl_mat[i,][which(read_counts[i,]!=-9)] <- as.vector(na.omit(pl[i,]))
}

print("formatted matrix")

# Replace NAs w/ the appropriate amount of zeros 
pl_mat_dat <- as.data.frame(pl_mat) # convert to df
pl_mat_dat[is.na(pl_mat_dat)] <- paste(as.character( rep(x = 0, times = p2+1) ), sep="' '", collapse=",")

print("replaced NAs w/ zeros")

# > Read in second PL file ----
pl2 <- readLines(pl_file2) %>% str_split(pattern = "\t", simplify = TRUE) 
pl2 <- apply(pl2, 2, function(x) gsub("^$|^ $", NA, x)) # replace blanks with NA

print("read in second file")

# load EBG input -- the total read counts matrix
read_counts_file2 <- read.csv(reads_file2, sep = "\t", header = F )
read_counts2 <- t(read_counts_file2)

# > Edit second file ----
# if values in row != -9,  assign PL to those spots in the new array [R, C] = [loci, individuals]
pl_mat2 <- matrix(data = NA,
                  nrow = dim(read_counts2)[1],
                  ncol = dim(read_counts2)[2] ) # make an empty matrix

l = dim(read_counts)[1]
for (i in 1:l){
  pl_mat2[i,][which(read_counts2[i,]!=-9)] <- as.vector(na.omit(pl2[i,]))
}

print("formatted second matrix")

# Replace NAs w/ the appropriate amount of zeros 
pl_mat2_dat <- as.data.frame(pl_mat2) # convert to df
pl_mat2_dat[is.na(pl_mat2)_dat)] <- paste(as.character( rep(x = 0, times = p2+1) ), sep="' '", collapse=",")

print("replaced NAs with zeros")

# > Add genotype information ----
colnames(pl_mat_dat) <- unlist(genos$V1)
colnames(pl_mat2_dat) <- unlist(genos2$V1)

# > bind the two dataframes together ----
#pl_all <- cbind(pl_mat, pl_mat2)
pl_all <- dplyr::bind_cols(pl_mat_dat, pl_mat2_dat)


# > Add SNP position information ----
rownames(pl_all) <- unlist(pos)

# > Export matrix ----
write.table(pl_all, file = paste0(fname, "-GL.txt"), quote = F)
