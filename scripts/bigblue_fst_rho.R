### Title: Estimating population differentiation 
### Author: Alyssa Phillips
### Date: 12/29/22

library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(stringr)
library(geosphere)

###
### Load data ----
###

# 1. Single-read genotypes from ANGSD ----
# g <- read.table(gzfile("~/Andropogon/pca/lowcov/all.andro.lowcov.50k.ibs.gz"), row.names = NULL, fill = T, header = T)
g <- read.table("~/Andropogon/kinship/all.andro.lowcov.100k.miss20.ibs", row.names = NULL, fill = T, header = T)

# Assign missing data to NA
g[g == -1] <- NA

# Transform so genotypes are columns
gt <- t(g[,5:161])

# 2. Metadata ----
meta <- read.csv("~/Andropogon/pca/lowcov/geno_meta_09212023.csv", header = T)
meta$realploidy <- as.factor(meta$realploidy)
meta$genotype <- as.factor(meta$genotype)

# 3. Make df ----
# Exclude clones
clones_to_drop = c("Afton_TX_3", "Afton_TX_1")
clones_to_drop_index = meta$index[meta$geno_short %in% clones_to_drop]

meta_noclones <- meta[!meta$geno_short %in% clones_to_drop,]
dim(meta_noclones)

gt <- gt[!rownames(gt) %in% clones_to_drop_index,]
dim(gt)

# Exclude 9x # IF E V W
meta_no9x <- meta_noclones[meta_noclones$ploidy == "2",]
meta_no9x <- meta_no9x[!is.na(meta_no9x$index),]
dim(meta_no9x)
gt <- gt[rownames(gt) %in% meta_no9x$index,]
dim(gt)

# Apply population and ploidy to single read genotypes 
east_pops <- c("TWE", "SUT", "FUL", "DES", "HAL", "MAN", "ESL", "CUI", "MON", "WEK", "WCG", "WIC", "KEN", "WAL" )
west_pops <- c("BOU","KON","MIL","AUS", "SAL","CDB","WEB","REL","BAR","AFT" )

ggroup <- lapply(meta_no9x$pop_short, function(x) ifelse(x %in% east_pops, "East", "West")) %>% unlist

df <- data.frame(
  # population = meta_noclones$population, # IF PAIRWISE
                 # pop_short = meta_noclones$pop_short, # IF PAIRWISE
                 # ploidy = meta_noclones$ploidy,# IF PAIRWISE
                 # pop.ploidy = paste0(meta_noclones$pop_short,".",meta_noclones$realploidy),# IF PAIRWISE
                 population = meta_no9x$population,  # IF E V W
                 pop_short = meta_no9x$pop_short,  # IF E V W
                 ploidy = meta_no9x$ploidy, # IF E V W
                 pop.ploidy = paste0(meta_no9x$pop_short,".",meta_no9x$realploidy), # IF E V W
                 ggroup = ggroup, # IF E V W
                 # genotype = meta$genotype,
                 gt
            )
df$pop_short <- as.factor(df$pop_short)
df$pop.ploidy <- as.factor(df$pop.ploidy)
dim(df)
df[1:5,1:5]


# ONLY IF Pairwise population comparison ----
# Drop populations with less than 3 genotypes per group
counts <- table(df$ploidy, df$pop_short)

tf_counts <- counts < 3

less_than_3 <- apply(tf_counts, 1, function(x) paste( names(which(x))) )
drop_pops <- c( paste0(less_than_3$'2', ".6"), paste0(less_than_3$'3', ".9") )

df <- dplyr::filter(df, !pop.ploidy %in% drop_pops)
df <- dplyr::filter(df, !ploidy %in% c(NA)) # ploidy only
dim(df) 

###
### Fst ----
###

# 1. Calculate allele frequencies for each sub-population and ploidy ----
# Controlling for missing data at each site
ps <- df %>% 
  group_by(ggroup) %>% # IF E v W
  # group_by(pop.ploidy) %>% # IF ALL PAIRWISE
  # slice_sample(n = 3) %>% # IF ALL PAIRWISE subsample to 3 genotypes per group
  group_map(~ colSums(.x[,5:dim(.x)[2]], na.rm = T) /  (dim(.x)[1] - colSums(is.na(.x[,5:dim(.x)[2]])) ) ) # this has the correct denominator 

subset_list <- df %>%
  group_by(ggroup) %>% # IF E v W
  dplyr::select(ggroup) %>% # IF E v W
  # group_by(pop.ploidy) %>% # IF ALL PAIRWISE
  # slice_sample(n = 3) %>% # IF ALL PAIRWISE
  # dplyr::select(pop.ploidy) %>% # IF ALL PAIRWISE
  unique()

# hist(ps[[1]])

# 2. Replace invariant sites with NA ----
for (i in 1:length(ps)){
  ps[[i]][ps[[i]] == 0] <- NA
  ps[[i]][ps[[i]] == 1] <- NA
}

# hist(ps[[1]])

# 3. Calculate Hsn ----
# Calculate for each site within each sub-pop (n = # of subpops)
Hsn <- lapply(ps, function(p){ 2 * p * (1 - p) })
# hist(Hsn[[1]])

# 4. Calculate HT ----
# I. Nest and pivot - create a one row df where each column is a list containing a single element
# group_names <- df$pop.ploidy %>% unique()
group_names <- unlist(subset_list)
ps_list <- tibble(group = group_names,
                  "vector" = ps)

# II. Make 2 way combos of the rows
ps_lists_comb <- tidyr::expand(ps_list, 
                        nesting(group, vector),
                        nesting(group2 = group, vector2 = vector))

# III. Filter redundancies
ps_lists_comb <- ps_lists_comb[!duplicated( t( apply(ps_lists_comb[,c(1,3)], 1, sort) ) ),]

ps_lists_red <- ps_lists_comb %>%
  filter(group != group2) %>%  # remove self-comparisions
  arrange(group, group2) %>% # sort rows
  mutate(vars = paste0(group, "-", group2)) # rearrange columns
  # select(contains("group"), everything())

# IV. Map functions across pairs
# calculate mean allele freq for each subpopulation pair at each site
pairs_mean_p <- ps_lists_red %>% 
  mutate(mean_p = map2(.x = vector, 
                       .y = vector2,
                       .f = ~{rowMeans(cbind(.x, .y))}))

# V. Calculate HT for each subpopulation pair at each site
# this properly accounts for missing data because the subpopulation allele frequencies have already taken this into account (I think)
pairs_HT <- pairs_mean_p %>%  
  mutate(HT = map(.x = mean_p,
                  .f = ~{2 * .x * (1 - .x)} 
                    ))

# 5. Calculate Fst ----
# ( HT - Hs ) / HT
hsn_list <- tibble(group = group_names ,
                  "vector" = Hsn)

# All pairs
# make an empty matrix where rows = pairs and cols = snps
Fst <- matrix(NA, nrow = dim(pairs_HT)[1], ncol = length(ps[[1]]) + 1 )
Fst[,1] <- pairs_HT$vars
Fst[1:5,1:5]
dim(Fst)

# i = "AFT.6-SAL.6"
for (i in pairs_HT$vars){
  temp_fst = as.numeric()
  ht = pairs_HT$HT[pairs_HT$vars == i] 
  
  sub_pop1 = unlist(strsplit(i, "-", fixed = T))[1]
  hsn1 = hsn_list$vector[hsn_list$group == sub_pop1]
  
  sub_pop2 = unlist(strsplit(i, "-", fixed = T))[2]
  hsn2 = hsn_list$vector[hsn_list$group == sub_pop2] 
  
  hs = (unlist(hsn1) + unlist(hsn2))/2 # average expected heterozygosity for subpopulations
  # temp_fst[which(unlist(ht) == 0)] = 0 # division by zero should be zero
  # temp_fst[which(unlist(ht) != 0)] = ( unlist(ht)[which(unlist(ht)!=0)] - hs[which(unlist(ht)!=0)] ) / unlist(ht)[which(unlist(ht)!=0)]
  temp_fst = (unlist(ht) - hs) / unlist(ht)
  Fst[,2:dim(Fst)[2]][Fst[,1] == i] <- as.numeric(temp_fst)
}

# 6. Calculate average Fst ----
## PAIRWISE
rownames(Fst) <- Fst[,1]
Fst <- Fst[,-1]
Fst_mat <- matrix(as.numeric(Fst),
       ncol = ncol(Fst))
rownames(Fst_mat) <- rownames(Fst)
# Fst_mat[rownames(Fst_mat) == "AFT.6-SAL.6",] %>% head(2000)

avg_fst <- rowMeans(Fst_mat, na.rm = T)
# hist(avg_fst)

# na_count <- rowSums(is.na(Fst_mat))
# plot(avg_fst, na_count)
# cor(avg_fst, na_count)

mean(avg_fst)
sd(avg_fst)

## E v W
# dim(Fst)
# Fst <- Fst[,-1] %>% as.numeric()
# mean(Fst, na.rm = T)
# sd(Fst, na.rm = T)




# 7. Plot avg Fst ----
avg_fst_df <- data.frame(fst = avg_fst,
                         pop1 = pairs_HT$group,
                         pop2 = pairs_HT$group2)


## W to E
# populationOrder <- c("BOU.3","AUS.3","KON.3","BOU.2", "MIL.2", "SAL.2", "Afton_TX.2","CDB.2", "BB.2","WEB.2","REL.2", "Austin_TX.2", "ESL.2", "MBTGPP.2","Konza.2","Cuivre_RiverSP_MO.2", "FUL.2", "DES.2", "12mi.2", "WI.2", "WAL.2", "Suther.2")
# populationOrder <- c("AUS.9","BOU.9","KON.9","AUS.6","BOU.6","KON.6", "MIL.6", "SAL.6", "AFT.6","CDB.6", "BAR.6","WEB.6","REL.6",  "ESL.6", "MAN.6","CUI.6", "FUL.6", "DES.6", "TWE.6", "KEN.6", "WAL.6", "SUT.6") # all

# W only
populationOrder_W <- c("AUS.9","BOU.9","KON.9","AUS.6","BOU.6","KON.6", "MIL.6", "SAL.6", "AFT.6","CDB.6", "BAR.6","WEB.6","REL.6") # all
avg_fst_df$pop1 <- factor(avg_fst_df$pop1, levels = populationOrder_W)
avg_fst_df$pop2 <- factor(avg_fst_df$pop2, levels = populationOrder_W)
avg_fst_df <- avg_fst_df[complete.cases(avg_fst_df),]

# avg_fst_df[avg_fst_df$pop1 == "BAR.6" & !avg_fst_df$pop2 %in% c("WEB.6", "REL.6"),]
avg_fst_df[avg_fst_df$pop1 == "BAR.6" & !avg_fst_df$pop2 %in% c("WEB.6", "REL.6"), c("pop1", "pop2")] <- avg_fst_df[avg_fst_df$pop1 == "BAR.6" & !avg_fst_df$pop2 %in% c("WEB.6", "REL.6"), c("pop2", "pop1")] 
avg_fst_df[avg_fst_df$pop1 == "CDB.6" & !avg_fst_df$pop2 %in% c("WEB.6", "REL.6", "BAR.6"),c("pop1", "pop2")] <- avg_fst_df[avg_fst_df$pop1 == "CDB.6" & !avg_fst_df$pop2 %in% c("WEB.6", "REL.6","BAR.6"), c("pop2", "pop1")] 
avg_fst_df[avg_fst_df$pop1 == "REL.6" & avg_fst_df$pop2 == c("SAL.6"),c("pop1", "pop2")] <- avg_fst_df[avg_fst_df$pop1 == "REL.6" & avg_fst_df$pop2 == c("SAL.6"), c("pop2", "pop1")] 
avg_fst_df[avg_fst_df$pop1 == "REL.6" & avg_fst_df$pop2 == c("WEB.6"),c("pop1", "pop2")] <- avg_fst_df[avg_fst_df$pop1 == "REL.6" & avg_fst_df$pop2 == c("WEB.6"), c("pop2", "pop1")] 
avg_fst_df[avg_fst_df$pop1 == "AUS.6" & avg_fst_df$pop2 %in% c("KON.9", "BOU.9","AUS.9"),c("pop1", "pop2")] <- avg_fst_df[avg_fst_df$pop1 == "AUS.6" & avg_fst_df$pop2 %in% c("KON.9", "BOU.9","AUS.9"), c("pop2", "pop1")] 
avg_fst_df[avg_fst_df$pop1 == "BOU.6" & avg_fst_df$pop2 %in% c("KON.9", "BOU.9"),c("pop1", "pop2")] <- avg_fst_df[avg_fst_df$pop1 == "BOU.6" & avg_fst_df$pop2 %in% c("KON.9", "BOU.9"), c("pop2", "pop1")] 
avg_fst_df[avg_fst_df$pop1 == "KON.6" & avg_fst_df$pop2 == c("KON.9"),c("pop1", "pop2")] <- avg_fst_df[avg_fst_df$pop1 == "KON.6" & avg_fst_df$pop2 == c("KON.9"), c("pop2", "pop1")] 

## de-duplicate
# avg_fst_df_sub <- avg_fst_df[!duplicated(avg_fst_df$fst),]
# avg_fst_df[upper.tri(avg_fst_df)]

# temp <- as.data.frame(t(apply(avg_fst_df_sub, 1, sort)))
# rownames(temp) <- rownames(avg_fst_df_sub)
# temp$fst <- as.numeric(temp$fst)
# 
# temp$pop1 <- factor(temp$pop1, levels = populationOrder)
# temp$pop2 <- factor(temp$pop2, levels = populationOrder)

avg_fst_df %>%
  ggplot(aes(x = pop1, y = pop2, fill = fst)) +
  geom_tile()  +
  scale_fill_distiller(palette = "BuGn", direction = 1, name = "Fst", limits = c(0,0.1)) +
  theme_ipsum() +
  # geom_text(aes(label = round(fst, digits = 2)), color = "white", size = 3) +
  theme(axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.key.size = unit(0.25, "in")) +
  scale_x_discrete(position = "top") 


ggsave(paste0("~/Andropogon/Figures/FST_WtoE.",Sys.Date(),".tiff" ), width = 10, height = 7, units = "in", )

## West only graph

# west <- c("Boulder_CO.9","Austin_TX.9","Konza.9","Boulder_CO.6", "NM.6", "SAL.6", "Afton_TX.6","CDB.6", "BB.6","WEB.6","REL.6", "Austin_TX.6", "ESL.6", "MBTGPP.6","Konza.6")
# west <- c("BOU.9","AUS.9","KON.9","BOU.6", "MIL.6", "SAL.6", "AFT.6","CDB.6", "BAR.6","WEB.6","REL.6", "AUS.6", "ESL.6", "MAN.6","KON.6")

# avg_fst_df_west <- avg_fst_df %>%
#   filter(pop1 %in% west) %>%
#   filter(pop2 %in% west)
# 
# avg_fst_df_west$pop2 <- factor(avg_fst_df_west$pop2, levels = west)
# avg_fst_df_west$pop1 <- factor(avg_fst_df_west$pop1, levels = west)
# 
# ggplot(avg_fst_df_west, aes(x = pop1, y = pop2, fill = fst)) +
#   geom_tile()  +
#   scale_fill_distiller(palette = "BuGn", direction = 1, name = "Fst"
#   ) +
#   theme_ipsum() +
#   # geom_text(aes(label = round(fst, digits = 2)), color = "white", size = 3) +
#   theme(axis.text.x = element_text(angle = 45, hjust= 1, size = 10),
#         axis.text.y = element_text(size = 10),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         panel.grid.major = element_blank())  

# library(ggridges)
# avg_fst_df %>%
#   filter(pop1 %in% west) %>%
#   filter(pop2 %in% west) %>%
#   ggplot(aes(x = fst)) +
#   geom_density() +
#   theme_bw() +
#   labs(x = "Fst", y = "Density") +
#   theme(axis.text = element_text(size = 10))

# 8. Test if 6x-9x and 9x-9x and 6x-6x Fst values are statistically different ----

avg_fst_df$ploidy1 <- stringr::str_split(avg_fst_df$pop1, pattern = "[.]", simplify = T)[,2]
avg_fst_df$ploidy2 <- stringr::str_split(avg_fst_df$pop2, pattern = "[.]", simplify = T)[,2]

avg_fst_df$test_type <- ifelse(test = avg_fst_df$ploidy1 == avg_fst_df$ploidy2,
       yes = ifelse(test = avg_fst_df$ploidy1 == "6",
                    yes = "hex",
                    no = "enn"),
       no = "mixed")
head(avg_fst_df)

## Compare means of 3 groups with Tukey test
boxplot(fst ~ test_type, data = avg_fst_df)

fst.mod <- lm(fst ~ test_type, data = avg_fst_df)

fst_means <- emmeans(fst.mod, specs = "test_type")
fst_diffs <- summary(contrast(fst_means, "pairwise"))
multcomp::cld(fst_means)

###
###  Rho ----
###

# Rho statistics developed Ronfort et al. (1998) 
# Explanation of calculation from Meirmans et al. (2018) Appendix 1

# 1. Calculate Hspn per pop ----
# Hspn = 1 - (p2 + q2)*((k-1) / k)
group_ploidies <- as.numeric(str_split(group_names, pattern = "[.]", simplify = T)[,2]) / 3 # IF PAIRWISE
group_ploidies <- c(2,2) # IF E V W
Hspn <- mapply(function(hsn,k){ hsn * (( k - 1) / k) }, 
               Hsn,
               group_ploidies,
               SIMPLIFY = F)
## allele frequencies of 0 will result in a Hspn of 0.5

# 2. Calculate Rho ----
## Within this calculation, calculate Hsp
# ( HT - Hs ) / (HT - Hsp)
hspn_list <- tibble(group = group_names ,
                   "vector" = Hspn)
# All pairs
# make an empty matrix where rows = pairs and cols = snps
rho <- matrix(NA, nrow = length(pairs_HT$vars), ncol = length(ps[[1]]) + 1 )
rho[,1] <- pairs_HT$vars
# rho[1:5,1:5]

# i = "12mi.6-Konza.9" # example
for (i in pairs_HT$vars){
  temp_rho = as.numeric()
  ht = unlist(pairs_HT$HT[pairs_HT$vars == i])
  
  sub_pop1 = unlist(strsplit(i, "-", fixed = T))[1]
  hsn1 = hsn_list$vector[hsn_list$group == sub_pop1]
  hspn1 = hspn_list$vector[hspn_list$group == sub_pop1]
  
  sub_pop2 = unlist(strsplit(i, "-", fixed = T))[2]
  hsn2 = hsn_list$vector[hsn_list$group == sub_pop2]
  hspn2 = hspn_list$vector[hspn_list$group == sub_pop2]
  
  hs = (unlist(hsn1) + unlist(hsn2))/2 # average expected heterozygosity for subpopulations
  hsp = (unlist(hspn1) + unlist(hspn2))/2 
  
  # temp_rho[which(ht==0)] = 0 # division by zero should be zero
  # temp_rho[which(ht!=0)] = ( unlist(ht)[which(unlist(ht)!=0)] - hs[which(unlist(ht)!=0)] ) / unlist(ht)[which(unlist(ht)!=0)]
  temp_rho = (ht - hs) / (ht - hsp)
  
  rho[,2:dim(rho)[2]][rho[,1] == i] <- as.numeric(temp_rho)
}

# 3. Calculate average Rho ----
## IF PAIRWISE
rownames(rho) <- rho[,1]
rho <- rho[,-1]
rho_mat <- matrix(as.numeric(rho),
                  ncol = ncol(rho))
rownames(rho_mat) <- rownames(rho)

avg_rho <- rowMeans(rho_mat, na.rm = T)


mean(avg_rho)
sd(avg_rho)

## E v W
dim(rho)
rho <- rho[,-1] %>% as.numeric()
mean(rho, na.rm = T)
sd(rho, na.rm = T)

# 4. Plot avg rho ---- 
avg_rho_df <- data.frame(rho = avg_rho,
                         pop1 = pairs_HT$group,
                         pop2 = pairs_HT$group2)

avg_rho_df$pop1 <- factor(avg_rho_df$pop1, levels = populationOrder_W)
avg_rho_df$pop2 <- factor(avg_rho_df$pop2, levels = populationOrder_W)
avg_rho_df <- avg_rho_df[complete.cases(avg_rho_df),]

avg_rho_df[avg_rho_df$pop2 == "REL.6" & !avg_rho_df$pop1 == "SAL.6",c("pop1", "pop2")] <- avg_rho_df[avg_rho_df$pop2 == "REL.6" & !avg_rho_df$pop1 == c("SAL.6"), c("pop2", "pop1")] 
avg_rho_df[avg_rho_df$pop2 == "BAR.6" & avg_rho_df$pop1 %in% c("AUS.9", "AUS.6"), c("pop1", "pop2")] <- avg_rho_df[avg_rho_df$pop2 == "BAR.6" & avg_rho_df$pop1 %in% c("AUS.9", "AUS.6"), c("pop2", "pop1")] 
avg_rho_df[avg_rho_df$pop2 == "WEB.6" & !avg_rho_df$pop1 %in% c("REL.6"),c("pop1", "pop2")] <- avg_rho_df[avg_rho_df$pop2 == "WEB.6" & !avg_rho_df$pop1 %in% c("REL.6"), c("pop2", "pop1")] 
avg_rho_df[avg_rho_df$pop2 == "CDB.6" & avg_rho_df$pop1 %in% c("AUS.9", "BOU.9","AUS.6","BOU.6"),c("pop1", "pop2")] <- avg_rho_df[avg_rho_df$pop2 == "CDB.6" & avg_rho_df$pop1 %in% c("AUS.9", "BOU.9","AUS.6","BOU.6"),c("pop2", "pop1")]
avg_rho_df[avg_rho_df$pop2 == "SAL.6" & !avg_rho_df$pop1 %in% c("CDB.6","BAR.6","WEB.6","REL.6"),c("pop1", "pop2")] <- avg_rho_df[avg_rho_df$pop2 == "SAL.6" & !avg_rho_df$pop1 %in% c("CDB.6","BAR.6","WEB.6","REL.6"),c("pop2", "pop1")]
avg_rho_df[avg_rho_df$pop2 == "MIL.6" & !avg_rho_df$pop1 %in% c("CDB.6","BAR.6","WEB.6","REL.6","SAL.6"),c("pop1", "pop2")] <- avg_rho_df[avg_rho_df$pop2 == "MIL.6" & !avg_rho_df$pop1 %in% c("CDB.6","BAR.6","WEB.6","REL.6","SAL.6"),c("pop2", "pop1")]
avg_rho_df[avg_rho_df$pop2 == "KON.6" & avg_rho_df$pop1 %in% c("AUS.9", "BOU.9","AUS.6","BOU.6"),c("pop1", "pop2")] <- avg_rho_df[avg_rho_df$pop2 == "KON.6" & avg_rho_df$pop1 %in% c("AUS.9", "BOU.9","AUS.6","BOU.6"),c("pop2", "pop1")]
avg_rho_df[avg_rho_df$pop2 == "BOU.6" & avg_rho_df$pop1 %in% c("AUS.9","AUS.6"),c("pop1", "pop2")] <- avg_rho_df[avg_rho_df$pop2 == "BOU.6" & avg_rho_df$pop1 %in% c("AUS.9","AUS.6"),c("pop2", "pop1")]
avg_rho_df[avg_rho_df$pop2 == "KON.9" & avg_rho_df$pop1 %in% c("AUS.9","BOU.9"),c("pop1", "pop2")] <- avg_rho_df[avg_rho_df$pop2 == "KON.9" & avg_rho_df$pop1 %in% c("AUS.9","BOU.9"),c("pop2", "pop1")]
avg_rho_df[avg_rho_df$pop2 == "BOU.9" & avg_rho_df$pop1 %in% c("AUS.9"),c("pop1", "pop2")] <- avg_rho_df[avg_rho_df$pop2 == "BOU.9" & avg_rho_df$pop1 %in% c("AUS.9"),c("pop2", "pop1")]



## All populations 
avg_rho_df %>%
  ggplot(aes(x = pop1, y = pop2, fill = rho)) +
  geom_tile()  +
  scale_fill_distiller(palette = "BuPu", direction = 1, name = "Rho", limits = c(0,0.15), labels = seq(0, 0.15, 0.05), breaks = seq(0, 0.15, 0.05)
  ) +
  theme_ipsum() +
  theme(axis.text.x = element_text(angle = 90, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.key.size = unit(0.25, "in"),
        panel.background = element_rect(fill = "transparent", colour = NA),  
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  scale_y_discrete(position = "right")


ggsave(paste0("~/Andropogon/Figures/Rho_WtoE.",Sys.Date(),".tiff" ), width = 10, height = 7, units = "in")

mean(avg_rho)
sd(avg_rho)

# 5. Test if 6x-9x and 9x-9x and 6x-6x Fst values are statistically different ----

avg_rho_df$ploidy1 <- stringr::str_split(avg_rho_df$pop1, pattern = "[.]", simplify = T)[,2]
avg_rho_df$ploidy2 <- stringr::str_split(avg_rho_df$pop2, pattern = "[.]", simplify = T)[,2]

avg_rho_df$test_type <- ifelse(test = avg_rho_df$ploidy1 == avg_rho_df$ploidy2,
                               yes = ifelse(test = avg_rho_df$ploidy1 == "6",
                                            yes = "hex",
                                            no = "enn"),
                               no = "mixed")
head(avg_rho_df)

## Compare means of 3 groups with Tukey test
boxplot(rho ~ test_type, data = avg_rho_df)

rho.mod <- lm(rho ~ test_type, data = avg_rho_df)

rho_means <- emmeans(rho.mod, specs = "test_type")
rho_diffs <- summary(contrast(rho_means, "pairwise"))
multcomp::cld(rho_means)


###
### (3) FST & RHO ----
###

# avg rho versus avg fst
plot(avg_rho, avg_fst)

# Histograms
hist(avg_rho,
     xlab = "Average Rho",
     main = NULL,
     xlim = c(0.04,0.1))
hist(avg_fst,
     xlab = "Average Fst",
     main = NULL,
     xlim = c(0.02,0.05))

###
### (4) Isolation by distance ----
###

# > Load population metadata ----
## For populations with two very close by groups, the average lat long value was used.
pop_meta <- read.csv("~/Andropogon/fst/all.andro.pop_meta.03032023.csv")
head(pop_meta)
str(pop_meta)
dim(pop_meta)

## Get correct population groupings
# avg_fst_df$popp1 <- str_split(avg_fst_df$pop1, "[.]", simplify = T)[,1] 
# avg_fst_df$popp2 <- str_split(avg_fst_df$pop2, "[.]", simplify = T)[,1] 
# head(avg_fst_df)

# > Calculate distance between points ----
pos <- as.matrix(cbind(pop_meta$Long, pop_meta$Lat))

dist <- distGeo(pos) # calculate distance in meters; distance on an ellipsoid
rownames(dist) <- unique(meta$pop_short)[which(pop_meta$population %in% unique(meta$population))]
colnames(dist) <- unique(meta$pop_short)[which(pop_meta$population %in% unique(meta$population))]
dist[1:5,1:5]

dist_km <- dist / 1000 # convert m to km

## Convert back into a data.frame
# dist_mat <- data.frame(pop1 = pop_meta$population[col(dist_km)], 
#                        pop2 = pop_meta$population[row(dist_km)], 
#                        distkm = c(dist_km),
#                        pair = paste0(pop_meta$population[col(dist_km)],"_",pop_meta$population[row(dist_km)])
#                        )
# head(dist_mat)

dist_unique <- which(upper.tri(dist_km, diag = TRUE), arr.ind = TRUE)
pop_names <- dimnames(dist_km)
dist_df <- data.frame(row = pop_names[[1]][dist_unique[, 1]],
           col = pop_names[[2]][dist_unique[, 2]],
           val = dist_km[dist_unique])
dist_df$p1_p2 <- paste0(dist_df$row, "-", dist_df$col)
dim(dist_df)

# > Plot FST vs Physical distance ----
avg_fst_df$fst_dst <- avg_fst_df$fst / (1 - avg_fst_df$fst)
# avg_fst_df$pairp <- paste0(avg_fst_df$pop1,"-",avg_fst_df$pop2) # pairs without ploidy
avg_fst_df$p1 <- stringr::str_split(avg_fst_df$pop1, pattern = "\\.", simplify = T)[,1]
avg_fst_df$p2 <-stringr::str_split(avg_fst_df$pop2, pattern = "\\.", simplify = T)[,1]
avg_fst_df$p1_p2 <- paste0(avg_fst_df$p1, "-", avg_fst_df$p2)

avg_fst_df_6 <- avg_fst_df[!grepl(pattern = "9", x = avg_fst_df$pairp),] # keep only 6x comparisons

head(dist_df)
head(avg_fst_df_6)
dist_df$p1_p2 %in% avg_fst_df_6$p1_p2 

### ok this whole thing needs to be re-thought - need to pull the pop pairs from avg_fst way early and then use it to select comparisions because order matters

unique(avg_fst_df$pop1) %>% length
dim(avg_fst_df)
## New matrix with fst distance and physical distance


dist_sorted <- dist_mat[match(avg_fst_df$pairp, dist_mat$pair),] %>% na.omit() # sort
ibd_df <- data.frame(fst = avg_fst_df$fst_dst, # combine
                dist = dist_sorted$distkm,
                pop1 = dist_sorted$pop1,
                pop2 = dist_sorted$pop2)
head(ibd_df)

## Plot
ggplot(as.data.frame(ibd_df), aes(x = dist, y = fst)) +
  geom_point() +
  xlab("Distance (m)") +
  ylab("Fst / (1 - Fst)")

cor(ibd_df$fst, ibd_df$dist) # Not correlated
