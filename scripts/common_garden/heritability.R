### Title: Heritability
### Author: Code written by Dan Runice and modified by Alyssa Phillips

devtools::install_github('deruncie/MegaLMM')
library(MegaLMM)
library(dplyr)
library(tidyverse)
library(ggplot2)

# (1) Load data ----
# All traits ----
## Created in linear_mod_all.R

## Leaf level traits
# df_lf <- read.csv("data/common_garden/leaf_traits.2024-03-27.csv", header = T)
df_lf <- read.csv("~/cg_andropogon/data/leaf_traits.2024-03-27.csv", header = T)
df_lf[,1:7] <- lapply(df_lf[,1:7], as.factor)
str(df_lf)

## Plant level traits
# df_pl <- read.csv("data/common_garden/plant_traits.2024-03-27.csv", header = T)
df_pl <- read.csv("~/cg_andropogon/data/plant_traits.2024-03-27.csv", header = T)
df_pl[,1:9] <- lapply(df_pl[,1:9], as.factor)
str(df_pl)

# kinship matrix ----
# k <- read.table("data/common_garden/all.andro.lowcov.100k.miss20.centereddiag.2023-09-27.txt")
k <- read.table("~/Andropogon/kinship/all.andro.lowcov.100k.miss20.centereddiag.2023-09-27.txt")
colnames(k) <- rownames(k)
k <- as.matrix(k)
k[1:5,1:5]
dim(k)

# (2) Prep data ----
# 1. Subset to only the 6x genotypes
# Subset k to 6x
hex_genos <- unique(df_pl$genotype[df_pl$ploidy == "6x"])

k_temp <- k[rownames(k) %in% hex_genos,]
k_hex <- k_temp[,colnames(k_temp) %in% hex_genos]
dim(k_hex)
rm(k_temp, k)

# 2. Subset df to genotypes in the kinship matrix
df_pl_sub <- filter(df_pl, genotype %in% rownames(k_hex))
df_pl_sub$genotype <- droplevels(df_pl_sub$genotype, exclude = "BL")
paste0("number of 6x genos: ", length(unique(df_pl_sub$genotype)))
rm(df_pl)

df_lf_sub <- filter(df_lf, genotype %in% rownames(k_hex))
df_lf_sub$genotype <- droplevels(df_lf_sub$genotype)
paste0("number of 6x genos: ", length(unique(df_lf_sub$genotype)))
rm(df_lf)

# (3) Make output data matrix ----
# rownames(blups) == sort(rownames(K)) # TRUE
# h2_mat <- matrix(NA, nrow = 15, ncol = 3) # 15 traits by 3 values (h2 + 95% CI)
# h2_df <- as.data.frame(h2_mat)
# rownames(h2_df) <- c("basalgrowthLOG", "relgrowth", "tillgrowth", "heightgrowth",
#                        "stomdenyr2", "stomlenyr2", "stomspiyr2", "ldmcyr2",
#                        "height", "tillers", "preretillersEXP",
#                        "slaINT", "length", "width", "thick")
# colnames(h2_df) <- c("h2", "lower_CI", "upper_CI")

## OR Reload output data matrix ----
h2_df <- read.csv("~/cg_andropogon/gmatrix/h2_estimates.2024-04-03.csv")
h2_df <- column_to_rownames(h2_df, var = 'X')
h2_df

# (4) Heritability: Run MegaLMM ----

## Single year traits ----
# df_hex = df_pl_sub # pick df
# single_year <- c("basalgrowthLOG", "relgrowth", "tillgrowth", "heightgrowth") # df_pl

# df_hex = df_lf_sub # pick df
# single_year <- c("stomspiyr2") # df_lf
# single_year <- c("stomdenyr2") # df_lf
# single_year <- c("stomlenyr2") # df_lf
# single_year <- c("ldmcyr2")

# Y = scale(df_hex[,single_year]) # traits only matrix
# run_ID = 'MegaLMM_run1_h2_singleyear_pl'


## Multi-year traits ----
# df_hex = df_pl_sub # pick df
# multi_year <- c("height")
# multi_year <- c("tillers")
# multi_year <- c("preretillersEXP") # df_pl

df_hex = df_lf_sub # pick df
multi_year <- c("slaINT")
# multi_year <- c("length")
# multi_year <- c("width")
# multi_year <- c("thick") # df_lf

### Reshape to add year column ----
scols <- c(paste0(multi_year, "yr1"), paste0(multi_year, "yr2"), "ploidy", "pop_short", "block", "genotype")
trait_wide <- df_hex %>%
  dplyr::select(all_of(scols)) %>%
  tidyr::pivot_longer(cols = c(paste0(multi_year, "yr1"), paste0(multi_year, "yr2")),
                      names_to = "year",
                      values_to = multi_year)
trait_wide <- trait_wide[complete.cases(trait_wide),]
trait_wide$year <- as.factor(trait_wide$year)

df_hex = trait_wide # re-assign name
rm(trait_wide)
head(df_hex)

Y = scale(df_hex[,multi_year]) # traits only matrix
run_ID = 'MegaLMM_run1_h2_multiyear'

# Specify parameters ----
run_parameters = MegaLMM_control(
  max_NA_groups = 3,
  scale_Y = FALSE,   # rescaled above
  h2_divisions = 20, # Each variance component is allowed to explain between 0% and 100% of the total variation. How many segments should the range [0,100) be divided into for each random effect?
  h2_step_size = NULL, # if NULL, all possible values of random effects are tried each iteration. If in (0,1), a new candidate set of random effect proportional variances is drawn uniformily with a range of this size
  burn = 100,  # number of burn in samples before saving posterior samples
  K = 10 # number of factors - might be way too big? Shouldn't matter much
)

rm(MegaLMM_state)
MegaLMM_state = setup_model_MegaLMM( ## CHANGE model type ----
  Y = Y,  
  # The n x p trait matrix
  # formula = ~ pop_short + (1|block) + (1|genotype),  # single year
  formula = ~ pop_short + (1|block) + (1|year) + (1|genotype),  # multi-year
  data = df_hex,         
  # the data.frame with information for constructing the model matrices
  relmat = list(genotype = k_hex), 
  # A list of covariance matrices to link to the random effects in formula.
  # each grouping variable in formula can be linked to a covariance matrix.
  # If so, every level of the grouping variable must be in the rownames of K.
  # additional rows of K not present in data will still be predicted 
  # (and therefore will use memory and computational time!)
  run_parameters = run_parameters,
  # This list of control parameters created above
  run_ID = run_ID
  # A run identifier. The function will create a folder with this name 
  # and store lots of useful data inside it
)
Lambda_prior = list(
  sampler = sample_Lambda_prec_horseshoe, 
  # function that implements the horseshoe-based Lambda prior 
  # described in Runcie et al 2020. 
  prop_0 = 0.1,    
  # prior guess at the number of non-zero loadings in the first and most important factor
  delta = list(shape = 3, scale = 1),    
  # parameters of the gamma distribution giving the expected change 
  # in proportion of non-zero loadings in each consecutive factor
  delta_iterations_factor = 100   
  # parameter that affects mixing of the MCMC sampler. This value is generally fine.
)
priors = MegaLMM_priors(
  tot_Y_var = list(V = 0.5,   nu = 5),      
  # Prior variance of trait residuals after accounting for fixed effects and factors
  # See MCMCglmm for meaning of V and nu
  tot_F_var = list(V = 18/20, nu = 20),     
  # Prior variance of factor traits. This is included to improve MCMC mixing, 
  # but can be turned off by setting nu very large
  h2_priors_resids_fun = function(h2s,n)  1,  
  # Function that returns the prior density for any value of the h2s vector 
  # (ie the vector of random effect proportional variances across all random effects. 
  # 1 means constant prior. 
  # n is the number of h2 divisions above (here=20)
  # 1-n*sum(h2s)/n linearly interpolates between 1 and 0, 
  # giving more weight to lower values
  h2_priors_factors_fun = function(h2s,n) 1, 
  # See above. 
  # sum(h2s) linearly interpolates between 0 and 1,
  # giving more weight to higher values
  # Another choice is one that gives 50% weight to h2==0: ifelse(h2s == 0,n,n/(n-1))
  Lambda_prior = Lambda_prior
  # from above
)
MegaLMM_state = set_priors_MegaLMM(MegaLMM_state, priors) # assign priors to model

# Manage missing data - drop as much as possible
maps = make_Missing_data_map(MegaLMM_state,max_NA_groups = ncol(Y)+1,verbose=F)
maps$map_results
MegaLMM_state = set_Missing_data_map(MegaLMM_state,maps$Missing_data_map_list[[1]]) ## CHANGE missing data selection ----

# Initialize the model object ----
estimate_memory_initialization_MegaLMM(MegaLMM_state) # estimate RAM required

MegaLMM_state = initialize_variables_MegaLMM(MegaLMM_state)
MegaLMM_state = initialize_MegaLMM(MegaLMM_state,verbose = T)

# Specify which parameters we want saved ----
## U = additive genetic values
## G = genetic covariances 
## R = residual covariances 
## h2 = narrow-sense (additive) heritability

MegaLMM_state$Posterior$posteriorSample_params
MegaLMM_state$Posterior$posteriorFunctions = list(
  # U = 'U_F %*% Lambda + U_R',
  # G = 't(Lambda) %*% diag(F_h2[1,]) %*% Lambda + diag(resid_h2[1,]/tot_Eta_prec[1,])',
  # R = 't(Lambda) %*% diag(1-F_h2[1,]) %*% Lambda + diag((1-resid_h2[1,])/tot_Eta_prec[1,])',
  h2 = '(colSums(F_h2[1,]*Lambda^2)+resid_h2[1,]/tot_Eta_prec[1,])/(colSums(Lambda^2)+1/tot_Eta_prec[1,])'
)
MegaLMM_state = clear_Posterior(MegaLMM_state) 

# Run Burnin period ----
n_iter = 100
for(i in 1:5) { 
  print(sprintf('Burnin run %d',i))
  # Factor order doesn't "mix" well in the MCMC.
  # We can help it by manually re-ordering from biggest to smallest
  MegaLMM_state = reorder_factors(MegaLMM_state,drop_cor_threshold = 0.6)
  # clear any previous collected samples because we've re-started the chain 
  MegaLMM_state = clear_Posterior(MegaLMM_state)
  # Draw n_iter new samples, storing the chain
  MegaLMM_state = sample_MegaLMM(MegaLMM_state,n_iter)
  # make diagnostic plots
  # traceplot_array(MegaLMM_state$Posterior$Lambda, name = file.path("~/cg_andropogon/gmatrix",MegaLMM_state$run_ID,'Lambda.pdf'))
  # traceplot_array(MegaLMM_state$Posterior$U,name =  file.path("~/cg_andropogon/gmatrix",MegaLMM_state$run_ID,'U.pdf'),
                  # facet_dim = 3)
  print(sprintf('Completed %d burnin samples', MegaLMM_state$current_state$nrun))
}
MegaLMM_state = clear_Posterior(MegaLMM_state)

# Collect posterior samples ----
## Use a thinning rate of 2, run for 1000 iterations -> 500 posterior samples
n_iter = 250
for(i in 1:4) {
  print(sprintf('Sampling run %d',i))
  MegaLMM_state = sample_MegaLMM(MegaLMM_state,n_iter) 
  MegaLMM_state = save_posterior_chunk(MegaLMM_state)
  print(MegaLMM_state)
}

# (4) Calculate predicted values and other posterior statistics ----
MegaLMM_state$Posterior = reload_Posterior(MegaLMM_state)

# Lambda_samples = load_posterior_param(MegaLMM_state,'Lambda')
MegaLMM_state$Posterior$posteriorSample_params = c(MegaLMM_state$Posterior$posteriorSample_params,'h2')
MegaLMM_state$Posterior$h2 = load_posterior_param(MegaLMM_state,'h2')

# Save model object
# saveRDS(MegaLMM_state,file = '~/cg_andropogon/gmatrix/MegaLMM_run1_h2_singleyear/MegaLMM_state_run_01_lf_singleyear.rds')

# Calculate posterior means
# Lambda = get_posterior_mean(MegaLMM_state,Lambda)
hcor = get_posterior_mean(MegaLMM_state, h2)
hcor

# Calculate highest posteriod density (HPD) intervals aka confidence intervals
## Default 95% CIs
h2_ci = get_posterior_HPDinterval(MegaLMM_state, h2)
h2_ci

# rownames(h2_df) %in% rownames(hcor)
# rownames(h2_df)[7] == rownames(hcor)


h2_df[multi_year,] <- data.frame(hcor = hcor[,1], ## CHANGE variable name ----
                          lower_CI = h2_ci[1,,],
                          upper_CI = h2_ci[2,,])
h2_df

# (5) Export estimates ----
write.csv(h2_df, paste0("~/cg_andropogon/gmatrix/h2_estimates.", Sys.Date(), ".csv"))
