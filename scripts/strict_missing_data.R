## Title: Strict missing data filter
## Author: Alyssa Phillips
## Date: 5/11/22

library("argparser", lib.loc = "R_libs")
library(stringr)

# > commandArgs ----
# give the arugment a name
ap <- arg_parser("grab SNPs with less than X percent missing data")

# add mandatory positional arguments (filename)
ap <- add_argument(ap, "gl", help = "GL matrix")

# parse arguments
argv <- parse_args(ap)

f <- as.character(argv$gl)
# f <- as.character("~/Andropogon/EBG/Chr05A-GL.txt")
fname <- strsplit(f, "-", )[[1]][1]

# Load GL matrix ----
gl_tab <- read.table(f)

# estimate the percent of individuals with no data per site
miss <- matrix(ncol=2, nrow=dim(gl_tab)[1])
colnames(miss) <- c("site", "per")

for (i in 1:dim(gl_tab)[1]){
  miss[i,2] <- sum(str_count(gl_tab[i,], '0,0,0|0,0,0,0')) / dim(gl_tab)[2] * 100
}

miss <- as.data.frame(miss)
miss$site <- rownames(gl_tab)

# plot missing data on the scaffold
#hist(miss$per)

# of snps under 5%
paste0("total number of snps before filtering:", dim(gl_tab)[1])
#paste0("number of snps under 5%:", sum(miss$per <= 5))
paste0("number of snps under 20%:", sum(miss$per <= 20))

# grab snps
#good <- gl_tab[miss$per <=5, ]
good <- gl_tab[miss$per <=20, ]
dim(good)

# export the good snps
#write.table(good, file = paste0(fname, ".miss5-GL.txt"), quote = F)
write.table(good, file = paste0(fname, ".miss20-GL.txt"), quote = F)
