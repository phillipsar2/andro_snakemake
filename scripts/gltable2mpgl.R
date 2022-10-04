### Title: PL table to mpgl format
### Author: Alyssa Phillips
### Date: 3/16/22

library(stringr)
library(MASS)
library("argparser", lib.loc = "R_libs")

ap <- arg_parser("Generate ENTROPY input files from GL matrix")

ap <- add_argument(ap, "--pl", help = "input PL table output from GATK")
ap <- add_argument(ap, "--out", help = "output prefix")

argv <- parse_args(ap)

# gl_file <- as.character("~/Andropogon/reports/filtering/all.AG.highcov.scaffold_32.PL.txt")
gl_file <- as.character(argv$pl)

# fname <- strsplit(gl_file, ".PL.txt")[[1]]
fname <- as.character(argv$out)

# (1) EBG PL to MPGL format -----------------------------------------------

print("Creating input genotype likelihood data for entropy...")

## >>  Read in GL file ----
gl_tab <- read.table(gl_file, header = T)

# >> Drop 9x sample ----
# 9x sample: INDE_PCRfree_Loretta_SAL-5_6_GCACGGAC_Andropogon_gerardii 
gl_sub <- subset(gl_tab, select = -INDE_PCRfree_Loretta_SAL.5_6_GCACGGAC_Andropogon_gerardii.PL)

# convert to matrix
#gl_mat <- unname(as.matrix(gl_tab[2:dim(gl_tab)[1], 2:dim(gl_tab)[2]]))
gl_mat <- unname(as.matrix(gl_sub[, 3:dim(gl_sub)[2]]))

# >> Split PLs by comma ----
# gl <- apply(gl_mat, 1:2, function(x){ as.numeric( unlist( strsplit(x, split = ",") ) ) })

ngls <- 3*dim(gl_mat)[2]
nsnp <- dim(gl_mat)[1]

gl <- matrix(NA, nrow = nsnp, ncol = ngls)

j = 1
k = 3
for (i in 1:dim(gl_mat)[2]){
  gl[ ,j:k] <- do.call(rbind, strsplit(gl_mat[ ,i], split = ",", fixed = TRUE))
  j = j + 3
  k = k + 3
}


# collapse 3-dim array to 2-dim matrix
# gl_collapse <- matrix(gl, ncol = prod(dim(gl)[1:2]), nrow = dim(gl)[3])
# dim(gl_collapse)
# gl_collapse[1:3,1:3]
# dim(gl) <- c(prod(dim(gl)[1:2]), dim(gl)[3])


# >> Prep metatdata ----
nind <- dim(gl_mat)[2]
nloci <- dim(gl_mat)[1]

genos <- colnames(gl_tab)[3:dim(gl_tab)[2]]
pos <- str_c(gl_tab$CHROM, "_", gl_tab$POS)
# pos <- rownames(gl_tab)



# >> Output mpgl ----
unlink(paste0(fname, ".mpgl")) # deletes the file with this name
cat(paste(nind, nloci, "\n"), file=paste0(fname, ".mpgl")) # print the nind and nloci into a file (line 1)
cat(paste((genos)), file=paste0(fname, ".mpgl"), append=T) # print genotype names in  line 2
for(locus in 1:nloci){
  cat(paste0("\n",pos[locus]," "), file=paste0(fname, ".mpgl"), append=T) # print locus name on line
  write.table(gl[locus,], file=paste0(fname, ".mpgl"), append=T, # write output table
              col.names=F, row.names=F, eol=" ", quote = F) # write gl values on same line
  # write.table(unlist(gl[locus,]), file=paste0(fname, ".mpgl"), append=T, # write output table
              # col.names=F, row.names=F, eol=" ", quote = F) # write gl values on same line
}

# >> Output ploidy input file ----
ploidy_mat <- apply(gl_mat, c(1,2), function(x){ length(unlist(str_split(x, pattern = ","))) - 1 })
tploidy_mat <- t(ploidy_mat) # flip to have individuals as rows and loci as columns

length(unlist(str_split(gl_mat[1,1], pattern = ",")))

write.table(ploidy_mat, file = paste0(fname, ".ploidy_inds.txt"), row.names=F, col.names=F, sep=" ")


# (2) MPGL to point estimates  -----------------------------------------------

# Create empty matrix for point estimates
mean_gl <- matrix(0, nrow = nloci, ncol = nind)

# >> Calculate point estimates of the genotypes ----
# multiplying genotypes (0,1,2) by the GLs = gives you a non-integer estimate of the genotype

# mean_gl <- apply(gl, c(1, 2), function(x){ temp <- sum(10^(-0.1 * unlist(x) )); sum((10^(-0.1 * unlist(x) )) * (0:mploidy)) / temp}) 

for (i in 1:nloci){
  for (j in 1:nind){
    g <- as.numeric(unlist(str_split(gl_mat[i,j], pattern = ",")))
    temp <- sum(10^(-0.1 * g))
    p <- length(g) - 1
    mean_gl[i,j] <- sum(( 10^(-0.1 * g) ) * (0:p)) / temp
  }
}

# >> Output point estimates ----
write.table(mean_gl, file=paste0(fname,"_pntest_meangl.txt"), row.names=F, col.names=F, quote=F, sep=" ")


# (3) Ancestry LDAK estimates  -----------------------------------------------

print("Obtaining admixture proportions for chain initialization...")

# PCA function
do.pca <- function(gmat){
  print("1. centering genotype matrix")
  gmn <- apply(gmat, 1, mean, na.rm=T) # mean at every site
  gmnmat <- matrix(gmn, nrow = nrow(gmat), ncol = ncol(gmat))
  gprime <- gmat - gmnmat # center around the mean
  
  print("2. finding pairwise covariance")
  gcovarmat <- matrix(NA, nrow = ncol(gmat), ncol = ncol(gmat))
  for(i in 1:ncol(gmat)){ # calculate covariance matrix
    for(j in i:ncol(gmat)){
      if (i == j){
        gcovarmat[i,j] <- cov(gprime[,i], gprime[,j], use="pairwise.complete.obs")
      }
      else{
        gcovarmat[i,j] <- cov(gprime[,i], gprime[,j], use="pairwise.complete.obs")
        gcovarmat[j,i] <- gcovarmat[i,j]
      }
    }
  }
  missing.inds <- unique( which( is.na(gcovarmat), arr.ind=T)[,2])
  print("3. final principal comp analysis")
  return(list(prcomp(x = na.omit(gcovarmat), center = TRUE, scale = FALSE), missing.inds))
}

# >> run PCA ----
return.val <- do.pca(mean_gl)

pcout <- return.val[1][[1]] # PCs
miss.inds <- return.val[2][[1]]
pcSummary <- summary(pcout) # summarize PCs

print("4. k-means, lda and writing initial admixture proportions to file")

# >> Estimate and output initial admixture prop ---- 

for(k in 2:14){
  init.admix <- matrix(0, nrow = nind, ncol=k)
  init.admix[miss.inds,] <- rep(1/k, k)
  
  kn <- kmeans(pcout$x[,1:5], k, iter.max = 10, nstart = 10, algorithm = "Hartigan-Wong") 
  ldakn <- lda(x = pcout$x[,1:5], grouping = kn$cluster, CV = T) 
  init.admix <- ldakn$posterior
  # output assignment probabilities to use as mean initialization values for q
  write.table(round(init.admix, 5), quote = F, row.names = F, col.names = F, file = paste0(fname,".qk",k,"inds.txt"))
}

