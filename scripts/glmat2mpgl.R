### Title: GL matrix to mpgl format
### Author: Alyssa Phillips
### Date: 3/16/22

library(stringr)
library(MASS, lib.loc = "R_libs")
library("argparser", lib.loc = "R_libs")

ap <- arg_parser("Generate ENTROPY input files from GL matrix")
 
ap <- add_argument(ap, "--gl", help = "input GL matrix")
ap <- add_argument(ap, "--out", help = "output prefix")
 
argv <- parse_args(ap)

#gl_file <- as.character("~/Andropogon/EBG/10k_lowcov-GL.txt")
gl_file <- as.character(argv$gl)

#fname <- strsplit(gl_file, "-GL.txt")[[1]]
fname <- as.character(argv$out)

# (1) EBG PL to MPGL format -----------------------------------------------

print("Creating input genotype likelihood data for entropy...")

## >>  Read in GL file ----
gl_tab <- read.csv(gl_file, sep = " ", header = T)

# convert to matrix
#gl_mat <- unname(as.matrix(gl_tab[2:dim(gl_tab)[1], 2:dim(gl_tab)[2]]))
gl_mat <- unname(as.matrix(gl_tab))

# >> Split PLs by comma ----
gl <- apply(gl_mat, 1:2, function(x){ as.numeric( unlist( strsplit(x, split = ",") ) ) })

# collapse 3-dim array to 2-dim matrix
# gl_collapse <- matrix(gl, ncol = prod(dim(gl)[1:2]), nrow = dim(gl)[2])
# dim(gl_collapse)

# >> Prep metatdata ----
nind <- dim(gl_mat)[2]
nloci <- dim(gl_mat)[1]

# genos <- gl_tab[1,1:dim(gl_mat)[2]]
genos <- colnames(gl_tab)
# pos <- gl_tab[2:dim(gl_mat)[1],1]
pos <- rownames(gl_tab)

# >> Output mpgl ----
unlink(paste0(fname, ".mpgl")) # deletes the file with this name
cat(paste(nind, nloci, "\n"), file=paste0(fname, ".mpgl")) # print the nind and nloci into a file (line 1)
cat(paste((genos)), file=paste0(fname, ".mpgl"), append=T) # print genotype names in  line 2
for(locus in 1:nloci){
  cat(paste0("\n",pos[locus]," "), file=paste0(fname, ".mpgl"), append=T) # print locus name on line
  # write.table(gl_collapse[locus,], file=paste0(fname, ".mpgl"), append=T, # write output table
  #             col.names=F, row.names=F, eol=" ") # write gl values on same line
  write.table(unlist(gl[locus,]), file=paste0(fname, ".mpgl"), append=T, # write output table
              col.names=F, row.names=F, eol=" ") # write gl values on same line
}

# >> Output ploidy input file ----
ploidy_mat <- apply(gl, c(1,2), function(x){ length( unlist(x) ) - 1 })
tploidy_mat <- t(ploidy_mat) # flip to have individuals as rows and loci as columns

write.table(ploidy_mat, file = paste0(fname, ".ploidy_inds.txt"), row.names=F, col.names=F, sep=" ")


# (2) MPGL to point estimates  -----------------------------------------------

# Create empty matrix for point estimates
mean_gl <- matrix(0, nrow = nloci, ncol = nind)

# >> Calculate point estimates of the genotypes ----
# multiplying genotypes (0,1,2) by the GLs = gives you a non-integer estimate of the genotype

# mean_gl <- apply(gl, c(1, 2), function(x){ temp <- sum(10^(-0.1 * unlist(x) )); sum((10^(-0.1 * unlist(x) )) * (0:mploidy)) / temp}) 

for (i in 1:nloci){
  for (j in 1:nind){
    temp <- sum(10^(-0.1 * unlist(gl[i,j])))
    p <- length(unlist(gl[i,j])) - 1
    mean_gl[i,j] <- sum(( 10^(-0.1 * unlist(gl[i,j])) ) * (0:p)) / temp
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

