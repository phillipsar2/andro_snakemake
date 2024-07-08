### Title: Transfer distance
### Author: Alyssa Phillips
### Date: 3/25/2024

library(dplyr)
library(ggplot2)
library(sommer)

# 0. Load data ----

## (1) MO data from field site 2021-2022 ----
# Extracted from http://agebb.missouri.edu/weather/history/report.asp?station_prefix=jef&start_month=5&end_month=9&start_day=1&end_day=30&start_year=2021&end_year=2022&period_type=1&convert=1&field_elements=70&field_elements=3&field_elements=51&field_elements=23&field_elements=48&field_elements=49&field_elements=73&field_elements=88&field_elements=146

col_env <- read.table("~/Andropogon/Common garden/columbia_env.txt", header = T)
col_env$MONTH <- as.factor(col_env$MONTH)
col_env$YEAR <- as.factor(col_env$YEAR)
head(col_env)
str(col_env)

## (2) Trait data ----
traits <- read.csv("~/cg_andropogon/data/plant_traits.2024-03-25.csv", header = T)
# head(traits)
# traits <- traits[!traits$plant == "BLANK",]
# traits$genotype[traits$genotype == "CAM1396"] <- "CAM_1396"
# traits$genotype <- as.factor(traits$genotype)

## (3) ClimateNA averages data ----
clm <- read.csv("~/cg_andropogon/scripts/Climate Data - ClimateNA/androclimatedata_popcodes.csv")
bou <- colMeans(clm[1:2,3:30]) # take the average environment from each boulder population
clm <- clm[2:16,] # drop one BOU row
clm[1,3:30] <- bou 

## (4) Kinship matrix ----
k <- read.table("~/Andropogon/kinship/all.andro.lowcov.100k.miss20.centereddiag.2023-09-27.txt")
colnames(k) <- rownames(k)
k <- as.matrix(k)
k[1:5,1:5]
dim(k)

## subset traits dataframe
traits_sub <- filter(traits, genotype %in% rownames(k)) %>%
  dplyr::select(relgrowth, genotype, block, ploidy, pop_short)
traits_sub$genotype <- as.factor(traits_sub$genotype)

traits_sub <- traits_sub[complete.cases(traits_sub),]
unique(traits_sub$genotype) %>% length()

traits_sub$genotype <- droplevels(traits_sub$genotype)
dim(traits_sub)

traits_sub[,2:5] <- lapply(traits_sub[,2:5], factor)

## subset kmat
K <- k[,colnames(k) %in% unique(traits_sub$genotype)]
K <- K[rownames(K) %in% unique(traits_sub$genotype),]
dim(K)
K[1:5,1:5]

## No clones need to be dropped from common garden dataset

## (5) Metadata ----
# meta <- read.csv("~/Andropogon/pca/lowcov/geno_meta_09212023.csv", header = T)
# head(meta)

## (6) Calculate BLUPS ----
rownames(K) %in% traits_sub$genotype
traits_sub$genotype %in%rownames(K) 
dim(K)

mod.blup <- mmer(relgrowth ~ 1 + ploidy,
                 random = ~vsr(genotype, Gu = K) + pop_short + block ,
                 rcov= ~units, 
                 nIters = 1000,
                 data = traits_sub, verbose = FALSE)

mod.blup$convergence # check convergence

grand_mean <- coef(mod.blup)[1,3] # intercept = grand mean?

blups <- as.data.frame(mod.blup$U$pop_short) + grand_mean # population BLUPs added to the grand mean
dim(blups)

# 1. Calculate climate variables for MO ----
## Mean annual temp (MAT) ----
colnames(col_env)

MAT <- (mean(col_env$AVG_TEMP_F) - 32) * (5/9) # Celcius
# aggregate(AVG_TEMP_F ~ YEAR, col_env, FUN = mean) * (5/9)

## Mean annual precipitation (MAP) ----
month_precip <- aggregate(TOTAL_PRECIP_IN ~ YEAR, col_env, FUN = sum)
MAP <- mean(month_precip$TOTAL_PRECIP_IN) * 25.4 # mm

month_temp <- aggregate(AVG_TEMP_F ~ MONTH + YEAR, col_env, FUN = mean)

# 2. Climate transfer distance ----
'We aim to test the association between transfer distance and mean relative growth estimated for each population using linear and quadratic regressions.' 

## (1) Estimate distance ----
env_df <- clm %>% 
  dplyr::select(ID2, MAT, MAP) 

env_df[env_df$ID2 =="Site",] <- c("Site", MAT, MAP) # replace with calculated values

env_df[,2:3] <- lapply(env_df[,2:3], as.numeric)

env_df$MAT_dist =   env_df[env_df$ID2 == "Site",]$MAT - env_df$MAT
env_df$MAP_dist =   env_df[env_df$ID2 == "Site",]$MAP - env_df$MAP


## (2) Model env dist ~ pop relgrowth ----
# e_df <- merge(traits_sub, env_df, by.x = "pop_short", by.y = "ID2")
# e_df

blups$pop_short <- stringr::str_split(rownames(blups), "pop_short", simplify = T)[,2]
e_df <- merge(blups, env_df[-16,], by.x = "pop_short", by.y = "ID2")
e_df

### MAT ----
mod.MAT <- lm(relgrowth ~ MAT_dist, e_df)
summary(mod.MAT)
anova(mod.MAT)
cor(e_df$relgrowth, e_df$MAT_dist)

### MAP ----
mod.MAP <- lm(relgrowth ~ MAP_dist, e_df)
summary(mod.MAP)
anova(mod.MAP)
cor(e_df$relgrowth, e_df$MAP_dist)


## (3) Plot predicted values ----

### MAT ----
p_mat <- sjPlot::plot_model(mod.MAT,
                   type = "pred", terms = "MAT_dist",
                   colors = c("black", "dark gray")
                   ) +
  theme_bw() +
  xlab(expression("MAT distance ("*degree*C*")")) +
  ylab("Mean change in aboveground biomass") +
  theme(plot.title = element_blank(),
        text = element_text(size = 12)) +
  geom_point(data = e_df, mapping = aes(x = MAT_dist, y = relgrowth), size = 2)

p_mat

ggsave(p_mat, file = paste0("~/cg_andropogon/Figures/MAT_transferdistance_scaled.", Sys.Date(),".jpeg"),
       width = 5,
       height = 5,
       units = "in")


### MAP ----
p_map <- sjPlot::plot_model(mod.MAP,
                   type = "pred", terms = "MAP_dist",
                   colors = c("black", "dark gray")) +
  theme_bw() +
  xlab("MAP distance (mm)") + ylab("Mean change in aboveground biomass") +
  theme(plot.title = element_blank(),
        text = element_text(size = 12)) +
  geom_point(data = e_df, mapping = aes(x = MAP_dist, y = relgrowth), size = 2)

p_map

ggsave(p_map, file = paste0("~/cg_andropogon/Figures/MAP_transferdistance_scaled.", Sys.Date(),".jpeg"),
       width = 4,
       height = 4,
       units = "in")

# 3. Geographic transfer distance ----
'We aim to test the association between transfer distance and mean relative growth estimated for each population using linear and quadratic regressions.' 

## (1) Estimate distance ----
# d_mat <- clm %>% 
#   dplyr::select(ID2, lat, long) %>%
#   stats::dist(method = "euclidean",
#               diag = TRUE, # print diagonal
#               upper = TRUE) %>% # print upper triangle
#   as.matrix()
# colnames(d_mat) <- clm$ID2
# rownames(d_mat) <- clm$ID2
# 
# d_site <- d_mat[-15,"Site"] # grab just the distance from the common garden 'Site'
# d_site <- d_site[order(names(d_site))]  %>% data.frame() # order
# d_site$pop_short <- rownames(d_site)
# colnames(d_site) <- c("distance", "pop_short")


## Haversine distance
library(geodist)
h_mat <- geodist(x = clm,measure = "haversine")
dim(h_mat)
colnames(h_mat) <- clm$ID2
rownames(h_mat) <- clm$ID2
h_site <- h_mat[1:14,"Site"]/1000  # kilometers
h_site 

blups$distance <- h_site[order(names(h_site))]
mod.geo <- lm(relgrowth ~ distance, blups)

blups

## (2) Model geo dist ~ pop relgrowth ----
# d_df <- merge(trait_means, d_site, by = "pop_short")
# d_df

# d_df <- merge(blups, d_site, by = "pop_short")
# d_df
# mod.geo <- lm(relgrowth ~ distance, d_df)

summary(mod.geo)
anova(mod.geo)
cor(blups$relgrowth, blups$distance)


## (3) Plot predicted values ----
# # predicted values - info located in the model summary
# plot(relgrowth ~ distance, data = d_df, type = "n", 
#      xlab = "Geographic distance", ylab = "Relative growth") # type n says don't add data, just axis
# points(relgrowth ~ distance, data = d_df, col = "black") 
# # could plot points seperately for each ploidy
# abline(a = coef(mod.geo)[1], b = coef(mod.geo)[2], col = "black" )

p_geo <- sjPlot::plot_model(mod.geo,
                   type = "pred",
                   colors = c("black", "dark gray"), 
                   terms = "distance"
                   ) +
  theme_bw() +
  xlab("Transfer distance (Km)") + 
  ylab("Mean change in aboveground biomass") +
  theme(plot.title = element_blank(),
        text = element_text(size = 12)) +
  geom_point(data = blups, mapping = aes(x = distance, y = relgrowth), size = 2) +
  scale_x_continuous(breaks = c(250,500,750,1000,1250,1500))

p_geo


ggsave(p_geo, file = paste0("~/cg_andropogon/Figures/geo_transferdistance_halversine_scaled.", Sys.Date(),".jpeg"),
width = 5,
height = 5,
units = "in")

## Combined figure for supplement
pp <- ggpubr::ggarrange(p_map, p_mat, ncol = 2,
                        labels = c("A", "B"),font.label = list(size = 22), hjust = -0.2)
ggsave(pp, file = paste0("~/cg_andropogon/Figures/map_and_mat_transferdist_scaled.", Sys.Date(),".jpeg"),
       width = 10,
       height = 5,
       units = "in")

# # 2. Env distance (PCs
# ## > PC space
# str(clim_var)
# 
# e_mat <- clim_var %>% 
#   dplyr::select(pop_short, PC1, PC3) %>%
#   unique() %>% 
#   stats::dist(method = "euclidean",
#               diag = TRUE, # print diagonal
#               upper = TRUE) %>% # print upper triangle
#   as.matrix()
# colnames(e_mat) <- clim_var$pop_short %>% unique()
# rownames(e_mat) <- clim_var$pop_short %>% unique()
# 
# e_site <- e_mat[-15,"Site"] # grab just the distance from the common garden 'Site'
# e_site <- e_site[order(names(e_site))] # order
# e_site
# 
# # > PC1
# rg <- aggregate(relgrowth ~ pop_short, df_pl, mean)
# str(clim_var)
# str(rg)
# pc_df <- merge(clim_var, rg, "pop_short")
# 
# # plot(x = pc_df$PC1, y = pc_df$relgrowth)
# 
# site_pc1 <- clim_var$PC1[clim_var$pop_short == 'Site']
# ggplot(pc_df, aes(x = PC1/site_pc1, y = relgrowth)) +
#   geom_smooth(method = "gam", se = T) +
#   geom_point() +
#   ylab("Population mean relative growth") + 
#   xlab("PC1 normalized by Site")
# 
# # > PC3
# 
# site_pc3 <- clim_var$PC3[clim_var$pop_short == 'Site']
# ggplot(pc_df, aes(x = PC3/site_pc3, y = relgrowth)) +
#   geom_smooth(method = "gam", se = T) +
#   geom_point() +
#   ylab("Population mean relative growth") + 
#   xlab("PC3 normalized by Site")
# 
# # 3. Plot
# rg$pop_short == names(e_site)
# 
# plot(y = rg$relgrowth, x = e_site)
# 
# td <- cbind(rg$relgrowth, 
#             e_site, 
#             d_site) %>%
#   as.data.frame()
# 
# mod.td <- lm(e_site ~ V1,td)
# summary(mod.td)
# car::Anova(mod.td)
# ggplot(td, aes(x = V1, y = e_site)) +
#   geom_smooth(method = lm, se = T) +
#   geom_point() +
#   xlab("Population mean relative growth") + ylab("Environmental distance from Columbia, MO (PC space)") +
#   # geom_text(hjust=0, vjust=2) +
#   theme_bw() +
#   annotate("text", label = "y ~ 5666 - 1649x", x = 3, y = 6000, size = 4) +
#   annotate("text", label = "p = 0.015", x = 3, y = 5500, size = 4)
# 
# mod.td <- lm(d_site ~ V1,td)
# summary(mod.td)
# car::Anova(mod.td)
# ggplot(td, aes(x = V1, y = d_site, label= rownames(td))) +
#   geom_smooth(method = lm, se = T) +
#   geom_point() +
#   xlab("Population mean relative growth") + ylab("Geographic distance from Columbia, MO") +
#   # geom_text(hjust=0, vjust=2) +
#   theme_bw() +
#   annotate("text", label = "y ~ 19.219 - 3.936x", x = 3, y = 20, size = 4) +
#   annotate("text", label = "p = 0.02082", x = 3, y = 19, size = 4)