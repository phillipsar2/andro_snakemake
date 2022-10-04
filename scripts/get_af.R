### Title: Calculate allele freqs from GL table
### Author: Alyssa Phillips, af funciton from Mitra Melon
### Date: 9/22/2022


library(stringr)

# Load data ----
pl <- read.table("~/Andropogon/ENTROPY/lowcov/10k.lowcov.miss5.fix-GL.txt", row.names = NULL, fill = T)

# remove duplicates
dups <- duplicated(pl$row.names)
pl_df <- pl[!duplicated(pl$row.names), ]

row.names(pl_df) <- pl_df$row.names

# drop row.names column
pl_df <- pl_df[1:dim(pl_df)[1], 2:dim(pl_df)[2]]

# Load metadata ----
meta <- read.table("~/Andropogon/ENTROPY/lowcov/sample_metadata.lowcov.txt")
colnames(meta) <- c("index", "geno", "pop", "ploidy")

# drop AN18N171_USPD16097227.GE07.AK12978_HJM33DSXX_L1 
meta <- meta[which(meta$geno != "AN18N171_USPD16097227.GE07.AK12978_HJM33DSXX_L1"),]

# (1) calculate allele frequencies from genotype likelihoods ----
# diploid_q = 1/2 * P(Aa) + P(aa)
# triploid_q = 1/3 * P(AAa) + 2/3 * P(Aaa) + P(aaa)

# write function
calcAF <- function(plTab){
  
  # This functional calculates the individual site allele frequency for diploids and triploids.
  # Code is based on function written by Mitra Melon
  # Input: PL table
  # Output: AF table
  
  # Convert PL to GP
  PL <- apply(plTab, c(1,2), function(x) 10 ^ ( as.numeric( str_split(x, ",", simplify = T) ) / -10 ) )
  gp <- apply(PL, c(1,2), function(x) unlist(x)/sum( unlist(x) ) )
  
  # Calculate allele freq
  af <- apply(gp, c(1,2), function(x){
    z <- unlist(x)
    if (length(z == 3)){ # diploids
      q <- 0.5 * z[2] + z[3]
      } else if (length(z)) { # triploid
        q <- 0.333 * z[2] + 0.666 * z[3] + z[4]
        } else {
        q <- NA # other; -9 will be the missing data symbol
        }
  }
  )
  
  return(af)
}

# Run function
bbAF <- calcAF(pl_df)

# Drop AN18N171_USPD16097227.GE07.AK12978_HJM33DSXX_L1
bbAF <- bbAF[,which(colnames(bbAF) != "AN18N171_USPD16097227.GE07.AK12978_HJM33DSXX_L1")]
dim(bbAF)
