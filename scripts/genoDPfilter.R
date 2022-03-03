# Title: Filtering SNPs by genotype depth
# Author: PoissionFilter function shared by Mitra Melon, script written by Alyssa Phillips
# Date: 3/2/2022

library(data.table)

# > Load vcf ----
# DF is the gdepth file outputted from vcftools
file <- as.character("reports/filtering/depth/lowcov/all.AG.lowcov.Chr07A.filtered.nocall.table")
P = 0.99
MISS = 0.2

#file <- as.character(snakemake@input[[1]])
# P = as.numeric(snakemake@input[[2]])
# MISS = as.numeric(snakemake@input[[3]])


fname <-  strsplit(file, split = ".table")

depth <- fread(file, header = T, sep = "\t") 
   
# # > Specify cutoff for qpois and % missingness
# p = 0.99
# miss = 0.2 # percent of missing data that is allowed

# > Filter data ----
PoissonFilter=function(DF,p,miss){
  
  #DF is the gdepth file outputted from gatk VariantsToTable with CHROM, POS, and FORMAT:DP fields
  #p is the cutoff for qpois
  #miss is the NA cutoff across all samples
  
  # prep depth data
  ORIG <- as.matrix(DF[ , -c(1:2)])
  info <- cbind( DF[ ,c(1:2)], chrPOS = paste(DF$CHROM,DF$POS, sep=":") )
  ORIG[ORIG == -1] <- NA
  
  # Calculate the depth cutoff for each genotype. The cutoff is the depth at which you can expect p percent of sites to have.
  cutoff <- apply(ORIG, 2, function(X) qpois(lambda = mean(X, na.rm=TRUE), p = p))
  
  # Mark sites NA if they exceed the maximum depth 
  for(c in 1:ncol(ORIG)){
    
    ID <- colnames(ORIG)[c] # ID of genotype
    
    out <- ORIG[ ,c] # depth vector for that genotype
    out[out > cutoff[ID]] <-  NA # assign NA to each site that has greater depth than the cutoff
    ORIG[ ,c] <- out # replace column with updated column **change to new matrix**
  }
  
  # Mark sites NA if the are less than the minimum depth of 1
  ORIG[ORIG < 1] <- NA
  
  # Count the number of NAs for each site
  avgNA = apply(ORIG, 1, function(X) sum( is.na(X) ))
  names(avgNA) <- info$chrPOS
  
  # Applying missingness filter
  D <- ncol(ORIG) * miss # calc the # of genos allowed to be missing
  keep <- avgNA[avgNA < D] # determine sites to keep
  keep <- info[info$chrPOS %in% names(keep), c(1:2)]
  
  # Print list of snps to keep
  return(keep)
  
  
}

good_snps <- PoissonFilter(DF = depth, p = P, miss = MISS)

# > Filtering results ----
paste0("Sites before filtering: ", dim(depth)[1])
paste0("Sites after filtering: ", dim(good_snps)[1])

# > Export table of good SNPs ----
write.table(good_snps, paste0(fname,".",P,"_", MISS,".txt"), sep = "\t", row.names = F, col.names = F, quote = F)
#write.table(good_snps, snakemake@output[[1]], sep = "\t", row.names = F, col.names = F, quote = F)

