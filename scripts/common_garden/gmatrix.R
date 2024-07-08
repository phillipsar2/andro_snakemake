### Title: G matrix
### Date: 7/5/23
### Author: Code written by Dan Runice and modified by Alyssa Phillips



# devtools::install_github('deruncie/MegaLMM')
library(MegaLMM)
library(sommer)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggnewscale)

# (1) Load data ----
# All traits ----
## Created in linear_mod_all.R
## Leaf-level traits are averages for the clone
df <- read.csv("~/cg_andropogon/data/all_traits.2024-03-26.csv", header = T)
df[,1:6] <- lapply(df[,1:6], as.factor)
str(df)

# kinship matrix ----
k <- read.table("~/Andropogon/kinship/all.andro.lowcov.100k.miss20.centereddiag.2023-09-27.txt")
colnames(k) <- rownames(k)
k <- as.matrix(k)
k[1:5,1:5]
dim(k)

# (2) Prep data ----
# 1. Subset to only the 6x genotypes
# Subset k to 6x
hex_genos <- unique(df$genotype[df$ploidy == "6x"])

k_temp <- k[rownames(k) %in% hex_genos,]
k_hex <- k_temp[,colnames(k_temp) %in% hex_genos]
dim(k_hex)

# 2. Subset df to genotypes in the kinship matrix
df_hex <- df[df$genotype %in% rownames(k_hex),]
dim(df_hex)
paste0("number of 6x genos: ", length(unique(df_hex$genotype)))

head(df_hex)
str(df_hex)

# 3. Make data matrix
Y = scale(df_hex[,-c(1:6)])
run_ID = 'MegaLMM_run1'

# (3) G matrix: Run MegaLMM ----

# Specify parameters ----
run_parameters = MegaLMM_control(
  max_NA_groups = 3,
  scale_Y = FALSE,   # should the columns of Y be re-scaled to have mean=0 and sd=1?
  h2_divisions = 20, # Each variance component is allowed to explain between 0% and 100% of the total variation. How many segments should the range [0,100) be divided into for each random effect?
  h2_step_size = NULL, # if NULL, all possible values of random effects are tried each iteration. If in (0,1), a new candidate set of random effect proportional variances is drawn uniformily with a range of this size
  burn = 100,  # number of burn in samples before saving posterior samples
  K = 20 # number of factors
)


MegaLMM_state = setup_model_MegaLMM(
  Y = Y,  
  # The n x p trait matrix
  formula = ~ pop_short + (1|block) + (1|genotype),  
  # We specify a fixed effect of population and a random effect for genotype and block
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
Lambda_prior = list( ## all default
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
MegaLMM_state = set_Missing_data_map(MegaLMM_state,maps$Missing_data_map_list[[6]])

# Initialize the model object ----
MegaLMM_state = initialize_variables_MegaLMM(MegaLMM_state)
MegaLMM_state = initialize_MegaLMM(MegaLMM_state,verbose = T)

# Specify which parameters we want saved ----
## U = additive genetic values
## G = genetic covariances 
## R = residual covariances 
## h2 = narrow-sense (additive) heritability

MegaLMM_state$Posterior$posteriorSample_params
MegaLMM_state$Posterior$posteriorFunctions = list(
  U = 'U_F %*% Lambda + U_R',
  G = 't(Lambda) %*% diag(F_h2[1,]) %*% Lambda + diag(resid_h2[1,]/tot_Eta_prec[1,])',
  R = 't(Lambda) %*% diag(1-F_h2[1,]) %*% Lambda + diag((1-resid_h2[1,])/tot_Eta_prec[1,])',
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
  traceplot_array(MegaLMM_state$Posterior$Lambda, name = file.path("~/cg_andropogon/gmatrix",MegaLMM_state$run_ID,'Lambda.pdf'))
  traceplot_array(MegaLMM_state$Posterior$U,name =  file.path("~/cg_andropogon/gmatrix",MegaLMM_state$run_ID,'U.pdf'),
                  facet_dim = 3)
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
G_samples = load_posterior_param(MegaLMM_state,'G')
R_samples = load_posterior_param(MegaLMM_state,'R')
dim(G_samples)
dim(R_samples)
Lambda_samples = load_posterior_param(MegaLMM_state,'Lambda')
MegaLMM_state$Posterior$posteriorSample_params = c(MegaLMM_state$Posterior$posteriorSample_params,'G','R','h2')
MegaLMM_state$Posterior$G = G_samples
MegaLMM_state$Posterior$R = R_samples
MegaLMM_state$Posterior$h2 = load_posterior_param(MegaLMM_state,'h2')

# Calculate posterior means
Lambda = get_posterior_mean(MegaLMM_state,Lambda)
Gcor = get_posterior_mean(MegaLMM_state, cov2cor(G)) # genetic correlation
Pcor = get_posterior_mean(MegaLMM_state, cov2cor(G+R)) # phenotypic correlation
hcor = get_posterior_mean(MegaLMM_state, h2)

# Calculate highest posteriod density (HPD) intervals aka confidence intervals
## Default 95% CIs
Gcors = get_posterior_HPDinterval(MegaLMM_state, cov2cor(G))
h2_ci = get_posterior_HPDinterval(MegaLMM_state, h2)




h2_df <- data.frame(trait = rownames(hcor),
                    hcor = hcor[,1],
           lower_CI = h2_ci[1,,],
           upper_CI = h2_ci[2,,])

# (5) Plot ----


# Plot G matrix ----
# Image(Gcor, dimnames=T, ) +
#   theme(panel.grid = element_blank())

Gcor2 <-
  Gcor %>%
  as.data.frame() %>%
  rownames_to_column("trait1") %>%
  pivot_longer(-trait1, names_to = "trait2", values_to = "value") %>%
  mutate(
    trait1 = factor(trait1),
    trait2 = factor(trait2)
  )

Gcors_lower <-
  Gcors[1,,] %>%
  as.data.frame() %>%
  rowid_to_column("trait1") %>%
  pivot_longer(-trait1, names_to = "trait2", values_to = "value") %>%
  mutate(
    trait1 = factor(trait1),
    trait2 = factor(gsub("V", "", trait2))
  )

Gcors_upper <-
  Gcors[2,,] %>%
  as.data.frame() %>%
  rowid_to_column("trait1") %>%
  pivot_longer(-trait1, names_to = "trait2", values_to = "value") %>%
  mutate(
    trait1 = factor(trait1),
    trait2 = factor(gsub("V", "", trait2))
  )

## Does the confidence interval cross zero?
Gcors_CI <- cbind(Gcors_lower$value, Gcors_upper$value)
cross_zero <- apply(cbind(Gcors_lower$value, Gcors_upper$value), 1, function(x) between(0, x[1], x[2] ))

Gcor2$cross_zero <- cross_zero

# Export G matrix
# write.csv(Gcor2, paste0("~/cg_andropogon/gmatrix/MegaLMM_run1/gmatrix_dataframe.", Sys.Date(), ".csv"), quote = F, row.names = F)

## Import G_df & h2 ----
Gcor2 <- read.csv("~/cg_andropogon/gmatrix/gmatrix_dataframe.2024-04-02.csv", header = T)
# h2 <- read.csv("~/cg_andropogon/gmatrix/h2_estimates.2024-04-03.csv", header = T)

# ordered and doubled for duplicate traits
h2_values <- c(0.33922451,0.33922451,0.12226719,0.28989426,0.28989426,0.22494398,0.13779780,0.24041588,0.19533394,0.19533394,0.05654885,0.05654885,0.16102520,0.16102520,0.04988974,0.04988974,0.08527293,0.08527293,0.18580503,0.14273880,0.11723638,0.27023863)

Gcor2$h2 <- c(NA)
Gcor2[Gcor2$trait1 == Gcor2$trait2,]$h2 <- h2_values
# Gcor2[Gcor2$trait1 == Gcor2$trait2,]$value <- h2_values # add heritabilities to diagnoal
Gcor2[Gcor2$trait1 == Gcor2$trait2,]$value <- NA
Gcor2[Gcor2$trait1 == Gcor2$trait2,]$cross_zero <- FALSE # all h2 signific

trait_order = c( "relgrowth", "heightgrowth", "tillgrowth", "basalgrowthLOG", "stomdenyr2", "stomlenyr2", "stomspiyr2", "ldmcyr2", "heightyr1", "tillersyr1", "preretillersEXPyr1", "lengthyr1", "widthyr1", "thickyr1", "slaINTyr1",  "heightyr2", "tillersyr2", "preretillersEXPyr2", "lengthyr2", "widthyr2", "thickyr2", "slaINTyr2" )
trait_labels= c( "RG", "HG", "TG", "BG", "SD", "SL", "SPI", "LDMC", "HY1", "TY1", "RTY1", "LY1", "WY1", "THY1", "SLAY1",  "HY2", "TY2", "RTY2", "LY2", "WY2", "THY2", "SLAY2" )



g_plot <- Gcor2 %>%
mutate(trait1 = factor(trait1, levels = trait_order),
       trait2 = factor(trait2, levels = trait_order)) %>%
ggplot(aes(trait1, trait2,
           # color = as.factor(cross_zero), 
           width=0.8, height=0.8), size = 3) +
  geom_tile(aes(fill = value)) +
  scale_fill_distiller(palette = "RdBu", limits=c(-1, 1),
                       na.value=NA, direction = 1) +
  labs(fill = bquote('r'['g'])) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        # legend.title = element_blank(),
        panel.grid = element_blank()) +
  # ggtitle("G matrix") +
  scale_colour_manual(name = "Significance", 
                      # breaks = names(df)[-1],
                      # labels = tempi_names,
                      values = c("black", "white")) +
  new_scale_fill() +
  geom_tile(aes(fill = h2), linewidth = NA) +
  scale_fill_distiller(palette = "YlGn" , na.value=NA, limits = c(0,0.4), direction = 1) +
  labs(fill = bquote('h'^2)) +
  geom_tile(data = Gcor2[Gcor2$cross_zero == FALSE,], 
            fill = NA, color = "black", linewidth = 0.5) +
  scale_x_discrete(labels = trait_labels) + 
  scale_y_discrete(labels = trait_labels)  
  # guides(fill= guide_legend(c("Gcor","h2")))
  # geom_segment(aes(x = 8.5, y = 8.5, xend = 8.5, yend = 22.5)) +
  # geom_segment(aes(x = 15.5, y = 8.5, xend = 15.5, yend = 22.5)) +
  # geom_segment(aes(y = 8.5, x = 8.5, yend = 8.5, xend = 22.5)) +
  # geom_segment(aes(y = 15.5, x = 8.5, yend = 15.5, xend = 22.5))
  

g_plot

ggsave(g_plot, paste0("~/cg_andropogon/Figures/Genetic_corr_plot",Sys.Date(),".jpeg"),
       height = 5,
       width = 6,
       units = "in")


## Plot a subset the matrix ----
trait_subset <- c("thickyr1", "widthyr1", "lengthyr1", "slaINTyr1")
trait_labels <- c("Thickness", "Width", "Length", "SLA")

Gcor_subset <- Gcor2 %>%
  mutate(trait1 = factor(trait1, levels = trait_order),
         trait2 = factor(trait2, levels = trait_order)) %>%
  dplyr::filter(trait1 %in% trait_subset & trait2 %in% trait_subset)
Gcor_subset$trait1 <- droplevels(Gcor_subset$trait1)
Gcor_subset$trait2 <- droplevels(Gcor_subset$trait2)
Gcor_subset$cross_zero[Gcor_subset$trait1 == Gcor_subset$trait2] <- TRUE

Gcor_subset %>%
  ggplot(aes(trait1, trait2,
             # color = as.factor(cross_zero), 
             width=0.9, height=0.9), size = 5) +
  geom_tile(aes(fill = value)) +
  scale_fill_distiller(palette = "RdBu", limits=c(-1, 1),
                       na.value=NA, direction = 1) +
  labs(fill = bquote('r'['g'])) +
  theme_minimal() +
  theme(text = element_text(size = 14),
    axis.text.x = element_text(angle = 90),
        axis.title = element_blank(),
        # legend.title = element_blank(),
        panel.grid = element_blank()) +
  # ggtitle("G matrix") +
  scale_colour_manual(name = "Significance", 
                      # breaks = names(df)[-1],
                      # labels = tempi_names,
                      values = c("black", "white")) +
  geom_tile(data = Gcor_subset[Gcor_subset$cross_zero == FALSE,], 
            fill = NA, color = "black", linewidth = 0.5) +
  scale_x_discrete(labels = trait_labels) +
  scale_y_discrete(labels = trait_labels)


# Image(hcor, dimnames=T)

# Plot h2 ----
ggplot(h2) +
  geom_segment( aes(x=trait, xend=trait, y=lower_CI, yend=upper_CI), color="grey") +
  geom_point( aes(x=trait, y=h2), color="black", size=2 ) +
  coord_flip()+
  theme_bw() +
  theme(
    legend.position = "none",
  ) +
  xlab("") +
  ylab(bquote('h'^2)) +
  scale_x_discrete(labels = rev(c("W", "TG","T","TH","SPI","SL","SD","SLA","RG","RT","L","LDMC","HG","H","BG"))) # flipping itself for some reason

ggsave(paste0("~/cg_andropogon/Figures/heritabilities.", Sys.Date(), ".jpeg"),
       width = 5,
       height = 6,
       unit = "in")

# Plot heatmap ----

## Load h2 ----
# h2 <- read.csv("~/cg_andropogon/data/trait_heritabilities_03282024.csv", header = F)

# demonominator of genetic correlation is the sqrt of heritabilities
# h2_combos <- data.frame(expand.grid(h2[,1], h2[,1]), expand.grid(h2[,2], h2[,2]))
# h2_combos$den <- h2_combos$Var1.1 * h2_combos$Var2.1
# ggplot(h2_combos, aes(x = Var1, y = sqrt(den), color= Var2)) +
#   geom_point() +
#   theme(axis.text.x = element_text(angle = 90))
# 
# out_df$gcor[out_df$Var1 == out_df$Var2] <- NA # set diagonal to NA

# trait_order = c( "relgrowth", "heightgrowth", "tillgrowth", "basalgrowthLOG", "stomdenyr2", "stomlenyr2", "stomspiyr2", "ldmcyr2", "heightyr1", "tillersyr1", "preretillersEXPyr1", "lengthyr1", "widthyr1", "thickyr1", "slaINTyr1",  "heightyr2", "tillersyr2", "preretillersEXPyr2", "lengthyr2", "widthyr2", "thickyr2", "slaINTyr2" )

# out_df %>%
# # mutate(Var1 = factor(Var1, levels = trait_order),
#        # Var2 = factor(Var2, levels = trait_order)) %>%
# ggplot(aes(Var1, Var2, fill= gcor)) + 
#   geom_tile() +
#   scale_fill_distiller(palette = "RdBu", limits=c(-1.6, 1)) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90),
#         axis.title = element_blank(),
#         legend.title = element_blank(),
#         panel.grid = element_blank())
# 
# ## Two traits against each other
# ggplot(df_hex, aes(x = heightgrowth, y = basalgrowthLOG, color = pop_short)) +
#   geom_point(alpha = 0.6) +
#      # xlab(paste0(trait1)) +
#      # ylab(paste0(trait2)) +
#   theme_bw() +
#   # ggtitle(paste0("rg = ", round(gcor,2))) +
#   geom_abline(intercept = 0, slope = 1)


# ggsave(paste0("~/cg_andropogon/data_explore/",trait1, "_", trait2, ".gcor.",Sys.Date(),".pdf"),
#        height = 6,
#        width = 8)