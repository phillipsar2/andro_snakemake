### Title: BEDASSLE analysis on big blue
### Author: Alyssa Phillips
### Date: 12/8/23

library(BEDASSLE)
library(dplyr)
library("argparser", lib.loc = "R_libs")

# input is a sorted vcf w/ header for a single ploidy level

# give the arugment a name
ap <- arg_parser("Running BEDASSLE v1")

# add mandatory positional arguments (filename)
ap <- add_argument(ap, "ac", help = "Allele count file; processed IBS file")

# additional arguments
ap <- add_argument(ap, "--gmeta", help = "Genotype meta data; one line per genotype")
ap <- add_argument(ap, "--pmeta", help = "Population meta data; one like per population")
ap <- add_argument(ap, "--env", help = "Environment data; one line per population")
ap <- add_argument(ap, "--pre", help = "Prefix for output files")
ap <- add_argument(ap, "--dir", help = "Output directory")


# parse arguments
argv <- parse_args(ap)

# Input
allele_counts <- as.character(argv$ac)
gmeta <- as.character(argv$gmeta)
pmeta <- as.character(argv$pmeta)
env_data <- as.character(argv$env)
pre <- as.character(argv$pre)
dir <- as.character(argv$dir)


# (0) Load data ----

## 1. Allele count data for all pops  ----
# Load all pops but test with a smaller amount first
ac <- read.table(allele_counts, row.names = NULL, fill = T, header = T)
ac_mat <- data.matrix(ac)

## 2. Load meta data ----
meta <- read.csv(gmeta, header = T)

colnames(ac_mat)[-c(1:4)] == meta$index
colnames(ac_mat)[-c(1:4)] <- meta$geno_short

### Exclude 9x ----
hexg <- meta$geno_short[meta$ploidy == 2]
hexg_pops <- meta$pop_short[meta$ploidy == 2]
ac_sub <- ac_mat[,colnames(ac_mat) %in% hexg]


### Remove invariant sites from subset (full dataset is ok) ----
ac_sub[ac_sub == -1 ] <- NA # set NAs

allele_count <- rowSums(ac_sub, na.rm = T)
paste0("fixed for ref allele: ", sum(allele_count == 0))

na_count <- rowSums(is.na(ac_sub)) # count NAs in dataset
paste0("fixed for alt allele: ", sum(allele_count == (dim(ac_sub)[2] - na_count) )) 

# remove sites fixed for alt and ref
allele_count[allele_count == 0] <- NA
allele_count[allele_count == (dim(ac_sub)[2] - na_count) ] <- NA
ac_sub_var <- ac_sub[!is.na(allele_count),]
paste0("Sites remaining: ", dim(ac_sub_var)[1]) # SNP count


### > Subset number of sites to 50K ----
ac_sub_sites <- ac_sub_var[1:50000,] # ALL IND

## 2. Geographic data ----
pop_meta <- read.csv(pmeta)

latlong <- pop_meta[,c(5,2:3)] %>% arrange(short_name)
colnames(latlong) <- c("pop", "lat", "long")

## 3. Environment data ----

### *Import env data ----
env_vars <- read.csv(env_data)

# (1) Calculate pairwise distances ----
# Standardize distances by dividing by SD

## 1. Distance matrix ----
D_mat <- stats::dist(latlong,
     method = "euclidean",
     diag = TRUE, # print diagonal
     upper = TRUE) # print upper triangle

D_mat <- D_mat/sd(D_mat)
D_mat <- as.matrix(D_mat)


## 2. Env matrix ----
E_mat_list <- list()
k = 1
for (j in env_vars_sub[,2:4]){
  e <- stats::dist(j,
              method = "euclidean",
              diag = TRUE, # print diagonal
              upper = TRUE) # print upper triangle
  e <- as.matrix(e)
  E_mat_list[[k]] <- e/sd(e)
  k = k + 1
}


# (2) Make input files ----
## 1. Allele count matrix ----
ac_t <- t(ac_sub_sites) %>% as.data.frame()
ac_t$pop <- hexg_pops[hexg %in% rownames(ac_t)]

counts <- aggregate(. ~ pop, data = ac_t, FUN = sum, na.rm = T,  na.action = NULL)

pop_names <- counts$pop # save pop names

counts <- counts %>%
  subset(select = -pop) %>% # drop pop names
  data.matrix() # convert to matrix
colnames(counts) <- NULL


## 2. Sample size matrix ----
# nrow = the number of populations
# ncol = the number of bi-allelic loci sampled
# Cells gives the number of chromosomes successfully genotyped at each locus in each population
pop_n <- data.frame(pop = ac_t$pop) %>% group_by(pop) %>% count() # pop sample sizes

sample_sizes <- matrix(data = NA, 
                       nrow = dim(counts)[1], 
                       ncol = (dim(counts)[2]))
paste0("Populations by number of sites: ", dim(sample_sizes)) # pops x sites

# Count the NAs from the popn at each site
count_nas <- aggregate(. ~ pop, data = ac_t, function(x) {sum(is.na(x))}, na.action = NULL)
count_nas <- subset(count_nas, select = -pop)

# Calculate final sample sizes
# (pop sample size * chromosomes sampled) - genotypes withoout data
for (site in 1:dim(sample_sizes)[2]){
  sample_sizes[,site] <- (as.numeric(pop_n$n) * 1) - count_nas[,site]
}


# (3) Build overdispersion model (MCMC_BB) ----

## 1. Set parameters ----
# * are set at default
delta  =  0.001 #*
aD_stp  =  0.0018 #*
aE_stp  =  0.04 #*
a2_stp  =  0.0035 #*
phi_stp  =  30 #*
thetas_stp =  0.07 #*
mu_stp  =  0.17 #*
ngen  =  1000000 # want closer to 10 million
printfreq  =  1000 # closer to 1000
savefreq  =  1000 # closer to 1000
samplefreq =  1000 # closer to 1000

#prefix  =  "Run1_ind3_10k_test"
#dir = "~/Andropogon/bedassle/"


nloci = dim(sample_sizes)[2]
npops = dim(sample_sizes)[1]

## 2. Run MCMC ---- 
MCMC_BB(counts = counts,
        sample_sizes = sample_sizes,
        D = D_mat,
        E = E_mat_list,
        k = npops,
        loci = nloci,
        delta = delta,
        aD_stp = aD_stp,
        aE_stp = aE_stp,
        a2_stp = a2_stp,
        thetas_stp = thetas_stp,
        phi_stp = phi_stp,
        mu_stp = mu_stp,
        ngen = ngen,
        printfreq = printfreq,
        savefreq = savefreq,
        samplefreq = samplefreq,
        directory = dir,
        prefix = pre)


