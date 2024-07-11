# Title: Evaluating effect of population and ploidy on plant functional traits
# Author: Alyssa Phillips

library(lme4)
library(lattice)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggpubr)
library(sommer)
library(emmeans)

###
### (0) Load data ----
###

# Load metadata ----
geno_meta <- read.csv("~/cg_andropogon/data/genotype_metadata.csv", header = T, skip = 0)
## add population shorthand
pop_short <- c("TWE", "AUS", NA, "KON", "SUT", "ESL", "MAN", "CDB", "FUL", "MIL", "BOU", "KEN", "REL", "SAL", "WAL")
pop_id <- data.frame("population" = unique(geno_meta$population), pop_short)
geno_meta <- merge(geno_meta, pop_id, "population")
head(geno_meta)

# Load trait data ----
## Leaf traits ----
leafyr1 <- read.csv("~/cg_andropogon/data/2021/2021 phenotype data - leaf traits_d04132022.csv", header = T)
leafyr1$sla <- (leafyr1$area_cm2 / 10000) / leafyr1$dry_mass_g # m2 / g
leafyr1$genoleaf <- paste0(leafyr1$genotype, "_", leafyr1$leaf)

lwyr2 <- read.csv("~/cg_andropogon/data/2022/2022 leaf length_width data - leaf traits_d03132023.csv", header = T)
lwyr2$genoleaf <- paste0(lwyr2$genotype, "_", lwyr2$leaf)

slayr2 <- read.csv("~/cg_andropogon/data/2022/2022 SLA_LDMC data - leaf traits_d03132023.csv", header = T)
slayr2$sla <- (slayr2$area_cm / 10000) / slayr2$dry_mass_g # m2/g
slayr2$ldmc <- (slayr2$dry_mass_g * 1000) / slayr2$fresh_mass_g

thickyr2 <- read.csv("~/cg_andropogon/data/2022/2022 phenotype data - leaf thickness_d03152023.csv", header = T)

## Basal area ----
bayr1 <- read.csv("~/cg_andropogon/data/2021/2021 phenotype data - basal aread03162022.csv", header = T)
bayr1$area <- bayr1$diameter_vertical * bayr1$diameter_horizontal * pi

bayr2 <- read.csv("~/cg_andropogon/data/2022/2022 phenotype data - basal area_d03132023.csv", header = T)
bayr2$area <- bayr2$diameter_vertical_cm * bayr2$diameter_horizontal_cm * pi

## Biomass ----
bioyr1 <- read.csv("~/cg_andropogon/data/2021/2021 phenotype data - biomass_d03102022.csv", header = F, skip = 2)
colnames(bioyr1) <- c("position", "block", "population","genotype", "biomass_g")
bioyr2 <- read.csv("~/cg_andropogon/data/2022/2022 phenotype data - biomass_d03132023.csv", header = T, skip = 0, nrows = 245)
colnames(bioyr2) <- c("position", "block", "genotype", "biomass_g")
bioyr2$bio_true <- as.numeric(bioyr2$biomass_g) - 61.875

relgrowth <- data.frame(relgrowth = bioyr2[order(bioyr2$position),]$bio_true / bioyr1[order(bioyr1$position),]$biomass_g, # yr2 / yr1
                        position = bioyr2[order(bioyr2$position),]$position)

## Plant height, survivorship, & tiller count ----
pstyr1 <- read.csv("~/cg_andropogon/data/2021/2021 phenotype data - survivorship & plant height & tillers_D03312022.csv", header = T)
pstyr2 <- read.csv("~/cg_andropogon/data/2022/2022 phenotype data - survivorship & plant height & tillers_d03132023.csv", header = F, skip = 2, nrows = 252)
colnames(pstyr2) <- c("position", "block","genotype", "survivorship", "plant_height_cm", "num_tillers", "num_reproductive_tillers")

## Stomata ----
# stomyr1 <- read.csv("~/cg_andropogon/data/2021/", header = T)

stomyr2 <- read.csv("~/cg_andropogon/data/2022/2022 stomata data - Stomata_d06162023.csv", header = T)
stomyr2[,6:12] <- lapply(stomyr2[,6:12], as.numeric)
stomyr2$meanlength <- rowMeans(stomyr2[,6:10], na.rm = T) # calculate mean length
stomyr2$density <- (stomyr2$stom_num/(stomyr2$area_um2/1000000)) # density = # stomata / area (mm^2)
stomyr2$spi <- stomyr2$density/sqrt(stomyr2$meanlength) # SPI = d/sqrt(mean length)

## 1. Whole plant traits dataframe ----
df_pl <- data.frame(position = pstyr1[order(pstyr1$position),]$position,
                    plant = pstyr1[order(pstyr1$position),]$genotype,
                    block = pstyr1[order(pstyr1$position),]$block,
                    genotype = pstyr1[order(pstyr1$position),]$genotype %>% stringr::str_sub(end = -4),
                    population = geno_meta[order(geno_meta$position),]$population,
                    pop_short = geno_meta[order(geno_meta$position),]$pop_short,
                    ploidy = geno_meta[order(geno_meta$position),]$ploidy,
                 survivoryr1 = pstyr1[order(pstyr1$position),]$survivorship,
                 survivoryr2 = pstyr2[order(pstyr2[,1]),]$survivorship,
                 heightyr1 = pstyr1[order(pstyr1$position),]$plant_height_cm,
                 heightyr2 = pstyr2[order(pstyr2$position),]$plant_height_cm,
                 heightgrowth = pstyr2[order(pstyr2$position),]$plant_height_cm / pstyr1[order(pstyr1$position),]$plant_height_cm,
                 tillersyr1 = pstyr1[order(pstyr1$position),]$num_tillers,
                 tillersyr2 = pstyr2[order(pstyr2$position),]$num_tillers,
                 tillgrowth = pstyr2[order(pstyr2$position),]$num_tillers /pstyr1[order(pstyr1$position),]$num_tillers, 
                 preretillersyr1 = pstyr1[order(pstyr1$position),]$num_reproductive_tillers/pstyr1[order(pstyr1$position),]$num_tillers,
                 preretillersyr2 = pstyr2[order(pstyr2$position),]$num_reproductive_tillers/pstyr2[order(pstyr2$position),]$num_tillers,
                 basalyr1 = bayr1[order(bayr1$position),]$area,
                 basalyr2 = bayr2[order(bayr2$position),]$area,
                 basalgrowth = bayr2[order(bayr2$position),]$area / bayr1[order(bayr1$position),]$area,
                 relgrowth = merge(pstyr1, relgrowth, "position", all.x = T, sort = T) %>% dplyr::select(relgrowth)
)

df_pl[,1:9] <- lapply(df_pl[,1:9], factor)
head(df_pl, 10)
str(df_pl)

### Remove outliers ----
df_pl$heightgrowth[df_pl$heightgrowth > 3] <- NA # height growth outliers
df_pl$basalgrowth[df_pl$basalgrowth > 25] <- NA # basal growth outliers
df_pl$basalgrowthLOG <- log(df_pl$basalgrowth)
df_pl$preretillersEXPyr1 <- (df_pl$preretillersyr1)^2
df_pl$preretillersEXPyr2 <- (df_pl$preretillersyr2)^2
df_pl$tillgrowth[df_pl$tillgrowth > 6] <- NA # tiller growth outliers

### Drop 2 poorly sequenced genotypes ----
#12MI-6_S2_L002_R2_001.fastq.gz
#CAM-1466_S7_L002_R2_001.fastq.gz

drop_g <- c("12mi_6", "CAM_1466")
df_pl <- df_pl[ !df_pl$genotype %in% drop_g ,]
dim(df_pl)

## 2. Leaf traits ----

pl <- c()
for ( i in 1:length(lwyr2$genoleaf)){
  pl[i] <- geno_meta[lwyr2$genotype[i] == geno_meta$plant, ]$ploidy
}

bl <- c()
for ( i in 1:length(lwyr2$genoleaf)){
  bl[i] <- geno_meta[lwyr2$genotype[i] == geno_meta$plant, ]$block
}


g <- c()
for ( i in 1:length(lwyr2$genoleaf)){
  g[i] <- geno_meta[lwyr2$genotype[i] == geno_meta$plant, ]$genotype
}

ps <- c()
for ( i in 1:length(lwyr2$genoleaf)){
  ps[i] <- geno_meta[lwyr2$genotype[i] == geno_meta$plant, ]$pop_short
}

df_lf <- data.frame(plant = lwyr2$genotype,
                    genotype = g,
                    leaf = lwyr2$leaf,
                    ploidy = pl,
                    block = bl,
                    pop_short = ps,
                    merge(lwyr2, leafyr1, "genoleaf", all.x = T, sort = F) %>% 
                      arrange(genoleaf) %>%
                      dplyr::select(population, length_cm.y, dry_mass_g, area_cm2, sla, leaf_thickness_mm, width_cm.y) %>% 
                      rename(lengthyr1 = length_cm.y,
                             dmassyr1= dry_mass_g,
                             slayr1 = sla,
                             thickyr1 = leaf_thickness_mm,
                             widthyr1 = width_cm.y,
                             areayr1 = area_cm2),
                    widthyr2 = lwyr2$width_cm,
                    lengthyr2 = lwyr2$length_cm,
                    dmassyr2 = slayr2[order(slayr2$genotype, slayr2$leaf),]$dry_mass_g,
                    fmassyr2 = slayr2[order(slayr2$genotype, slayr2$leaf),]$fresh_mass_g,
                    areayr2 = slayr2[order(slayr2$genotype, slayr2$leaf),]$area_cm,
                    slayr2 = slayr2[order(slayr2$genotype, slayr2$leaf),]$sla,
                    ldmcyr2 = slayr2[order(slayr2$genotype, slayr2$leaf),]$ldmc,
                    thickyr2 = thickyr2[order(thickyr2$genotype, thickyr2$leaf),]$leaf_thickness_mm,
                    denyr2 = (thickyr2[order(thickyr2$genotype, thickyr2$leaf),]$leaf_thickness_mm * 1000)^-1 * (slayr2[order(slayr2$genotype, slayr2$leaf),]$sla * 1000)^-1 * (slayr2[order(slayr2$genotype, slayr2$leaf),]$ldmc )^-1,
                    stomdenyr2 = stomyr2[order(slayr2$genotype, stomyr2$genotype),]$density,
                    stomlenyr2 = stomyr2[order(slayr2$genotype, stomyr2$genotype),]$meanlength,
                    stomspiyr2 = stomyr2[order(slayr2$genotype, stomyr2$genotype),]$spi
                    )

head(df_lf)
df_lf[,1:7] <- lapply(df_lf[,1:7], factor)
str(df_lf)

### remove outliers ----
df_lf$denyr2[df_lf$denyr2 > 2.5] <- NA
df_lf$stomdenyr2[df_lf$stomdenyr2 > 1500] <- NA
df_lf$stomspiyr2[df_lf$stomspiyr2 > 140] <- NA
df_lf$ldmcyr2[df_lf$ldmcyr2 > 1000] <- NA
df_lf$lengthyr2[df_lf$lengthyr2 > 30] <- NA

invNormTransform = function(x) {
  qnorm((rank(x,na.last="keep")-3/8)/sum(!is.na(x)))
}
df_lf$slaINTyr1 <- invNormTransform(df_lf$slayr1) # transform SLA
df_lf$slaINTyr2 <- invNormTransform(df_lf$slayr2)

### Drop 2 poorly sequenced genotypes ----
df_lf <- df_lf[ !df_lf$genotype %in% drop_g ,]
dim(df_lf)

# 3. Export  traits ----

# All traits - one value per plant - for genetic covariance estimates
## Average leaf-level traits across leaves
# trait_names <- colnames(df_lf[,8:25])
# df_lf_sum <- aggregate(cbind(lengthyr1, lengthyr2, slaINTyr1, slaINTyr2, thickyr1, thickyr2, widthyr1, widthyr2,  ldmcyr2, stomdenyr2, stomlenyr2, stomspiyr2)  ~ pop_short + plant + ploidy + block, df_lf, FUN = mean)
# 
# ## Subset df_pl to traits we care about
# df_pl_traits <- df_pl %>% dplyr::select(position, plant, block, genotype, pop_short, ploidy, heightyr1, heightyr2, heightgrowth, tillersyr1, tillersyr2, tillgrowth, basalgrowthLOG, relgrowth, preretillersEXPyr1, preretillersEXPyr2)

## Combine with whole plant trait data
# df_all <- merge(df_pl_traits, df_lf_sum, "plant", all.x = T, sort = F) %>%
#   dplyr::select(!c(pop_short.y, ploidy.y, block.y)) %>%
#   arrange(position) %>%
#   filter(!plant == "BLANK") # exclude blanks
# colnames(df_all)[c(3,5:6)] <- c("block", "pop_short","ploidy")
# head(df_all)

## Export
# write.csv(df_all, file = paste0("~/cg_andropogon/data/all_traits.",Sys.Date(),".csv"), quote = F, col.names = T, row.names = F)

### Leaf traits
# export
# write.csv(df_lf, file = paste0("~/cg_andropogon/data/leaf_traits.",Sys.Date(),".csv"), quote = F, col.names = T, row.names = F)

### Plant traits
# write.csv(df_pl, file = paste0("~/cg_andropogon/data/plant_traits.",Sys.Date(),".csv"), quote = F, col.names = T, row.names = F)

# 4. Load kinship matrix ----
# Centered diagonal - diagonal doesn't effect estimates
k <- read.table("~/Andropogon/kinship/all.andro.lowcov.100k.miss20.centereddiag.2023-09-27.txt")
colnames(k) <- rownames(k)
k <- as.matrix(k)
k[1:5,1:5]
dim(k)
## subset kmat
K <- k[,colnames(k) %in% unique(df_lf$genotype)]
K <- K[rownames(K) %in% unique(df_lf$genotype),]
dim(K)
K[1:5,1:5]

###
### (1) Check out your data ----
###

# Function to make exploratory plots for each variable
plot_data <- function(df, y_var, ord){
  trait_means <- aggregate(as.formula(paste0(y_var, "~ pop_short + plant + ploidy + block")), df, FUN = mean)
  genotype_means <- aggregate(as.formula(paste0(y_var, "~ pop_short + genotype + ploidy")), df, FUN = mean)
  
  p1 <- ggplot(genotype_means, aes(x = forcats::fct_reorder(.data[[ord]], .data[[y_var]], .fun = mean), y = .data[[y_var]], fill = ploidy)) + 
    geom_boxplot()+
    geom_point(position = position_jitterdodge(jitter.width = 0)) +
    xlab("Population") +
    ylab(paste0(y_var)) +
    theme_bw() +
    scale_fill_manual(values = c("white", "darkgray"))
  
  p2 <- ggplot(trait_means, aes_string(x = y_var)) +
    geom_histogram() + 
    theme_bw()
  
  p3 <- ggplot(trait_means, aes_string(x = "block", y = y_var)) +
    geom_boxplot() +
    geom_jitter(width = 0.2) +
    theme_bw()
  
  p4 <- genotype_means %>% 
    na.omit() %>%
    ggplot( aes_string(x = "ploidy", y = y_var)) +
    geom_boxplot() +
    geom_jitter(width = 0.2) +
    theme_bw()
  
  ggarrange(p1, 
            ggarrange(p2, p3, p4, ncol = 3),
            nrow = 2,
            heights = c(2.25, 2)
  )
  
  ggsave(paste0("~/cg_andropogon/data_explore/",y_var, ".explore.",Sys.Date(),".pdf"),
         height = 6,
         width = 8)
}

## 1. Plot data ----

## plots for whole-plant traits
# for (i in colnames(df_pl)[10:22]){
#   plot_data(df_pl, i, "pop_short")
# }

# plot_data(df_pl, "preretillersEXPyr2", "pop_short")

## plots for leaf traits
# for (i in colnames(df_lf)[8:21]){
#   plot_data(df_lf, i, "pop_short")
# }

# plot_data(df_lf, "slaINTyr2", 'pop_short')

## 2. Effect of year ----
## c(height, perretillers, length, width)
df_pl %>%
  dplyr::select(heightyr1, heightyr2, ploidy, pop_short, genotype) %>%
  tidyr::pivot_longer(cols = c(heightyr1, heightyr2),
                      names_to = "year",
                      values_to = "height") %>%
  aggregate(height ~ pop_short + genotype + ploidy + year, FUN = mean) %>% 
  na.omit() %>% 
  ggplot(aes(x = forcats::fct_reorder(pop_short, height), y = height, fill = year)) +
    geom_boxplot() +
    geom_point(position = position_jitterdodge(), size = 1) +
    xlab("Population") +
  ylab("Height (cm)") +
    theme_bw() +
  scale_fill_manual(values = c("white","lightblue"),
                    labels = c("1", "2"), 
                    name = "Year")

# 3.  Assess differences due to data collector ----
## Width 
aggregate(width_cm ~ genotype + name, lwyr2, FUN = mean) %>%
ggplot(aes(x = name, y = width_cm)) +
  geom_boxplot() +
  geom_jitter(width = .2)
t.test(lwyr2[lwyr2$name == "Nathan",4], lwyr2[lwyr2$name == "Alyssa",4])

## Length 
aggregate(length_cm ~ genotype + name, lwyr2, FUN = mean) %>%
  ggplot(aes(x = name, y = length_cm)) +
  geom_boxplot() +
  geom_jitter(width = .2)
t.test(lwyr2[lwyr2$name == "Nathan",5], lwyr2[lwyr2$name == "Alyssa",5])

## Dry mass
aggregate(dry_mass_g ~ genotype + namefordrymass, slayr2, FUN = mean) %>%
  ggplot(aes(x = namefordrymass, y = dry_mass_g)) +
  geom_boxplot() +
  geom_jitter(width = .2)
t.test(slayr2[slayr2$namefordrymass == "Nathan",5], slayr2[slayr2$namefordrymass == "Muyun",5])

## 4. Figures for publication ----

## * Stomata figures ----
genotype_means <- aggregate(stomdenyr2 ~ pop_short + genotype + ploidy, df_lf, FUN = mean) 

summary(genotype_means$ploidy) # sample size for plot

# y_name <- expression(paste("Mean stomata length (", mu, "m)"))


## Ploidy boxplot
axis_names <- c( expression(paste("6", italic("x"))),
                 expression(paste("9", italic("x")))
)
stom_p <- genotype_means %>% 
  na.omit() %>%
  ggplot(aes(x = ploidy, y = stomdenyr2, 
             # fill = ploidy
             )) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 2) +
  theme_bw() +
  xlab("Cytotype") +
  ylab(expression('Stomata density (m'^2*')')) +
  scale_x_discrete(labels = axis_names) +
  # scale_fill_manual(values = c("white", "#666666")) +
  guides(fill = "none") +
  theme(text = element_text(size = 12))

stom_p

ggsave(stom_p,
       file = paste0("~/cg_andropogon/Figures/stomatadensity_by_ploidy.",Sys.Date(),".jpeg"),
       width = 5,
       height = 5,
       units = "in")


# * Relative growth fig ----
genotype_means <- aggregate(relgrowth ~ pop_short + genotype + ploidy, df_pl, FUN = mean) 

# genotype_means[genotype_means$pop_short %in% c("BOU", "AUS", "KON"),] %>% 
genotype_means %>%
  ggplot(aes(x = ploidy, y = relgrowth, fill = ploidy)) +
  geom_boxplot() +
  geom_point(position=position_jitterdodge(jitter.width = 0.3), size = 1) +
  theme_bw() +
  xlab("Population") +
  ylab("Relative aboveground growth") +
  scale_fill_manual(values = c("white", "#666666")) +
  guides(fill = "none") +
  theme(text = element_text(size = 12))

# ggsave(
#        file = paste0("~/cg_andropogon/Figures/relgrowth_by_ploidy.",Sys.Date(),".jpeg"),
#        width = 4,
#        height = 4,
#        units = "in")

###  Trait v trait
trait1_means <- aggregate(stomdenyr2 ~ pop_short + genotype + ploidy + MAP, df_lf_clim, FUN = mean)
# biomass_means <- merge(bioyr2,geno_meta, by.x = "genotype", by.y = "plant") %>% 
#   aggregate(biomass_g ~ genotype.y, FUN = mean) 
# stom_means <- aggregate(stomlenyr2 ~ pop_short + genotype + ploidy, df_lf, FUN = mean)
trait2_means <- aggregate(relgrowth ~ pop_short + genotype + ploidy + MAP, df_pl_clim, FUN = mean)

# growth_stom_means <- merge(growth_means, stom_means, by = "genotype")
# growth_stom_means <- merge(stom_means, biomass_means, by.x = "genotype", by.y = "genotype.y")
combo_means <- merge(trait1_means, trait2_means, by = "genotype")

# ggplot(growth_stom_means, aes(x = relgrowth, y = stomlenyr2)) +
#   geom_point() 
#   # geom_smooth()
hex_colors <- c("#FF0000", "#00FF00", "#0066FF", "#FFFF00", "#00FFFF", "#FF00FF", "#FFA500", "#800080", "#008000", "#FF6347", "#00CED1", "#9A2FE2", "#7FFFD4", "#FFD700")

combo_means$group <- ifelse(combo_means$MAP.x > 750, "Wet", "Dry")

combo_means %>%
ggplot(aes(x = relgrowth, y = stomdenyr2, group = group)) +
  geom_point(aes(col = pop_short.x)) +
  geom_smooth(method = lm,se = F, color = "black") +
  theme_bw() +
  xlab("relgrowth") + ylab("relgrowth") +
  theme(text = element_text(size = 12)) +
  scale_color_manual(values = hex_colors, name = "Population")

## Population boxplot
# pop_p <- genotype_means %>% 
#   na.omit() %>%
#   ggplot(aes(x = forcats::fct_reorder(.data[['pop_short']], .data[['stomlenyr2']], .fun = mean), y = stomlenyr2, fill = ploidy)) +
#   geom_boxplot() +
#   geom_point(aes(fill = ploidy), size = 1, position = position_jitterdodge(), shape = 21) +
#   # geom_jitter(width = 0.1, size = 1, shape = 21) +
#   theme_bw() +
#   xlab("Population") +
#   ylab(y_name) +
#   theme(text = element_text(size = 12)) +
#   scale_fill_manual(values = c("white", "#666666")) +
#   guides(fill = "none")
# 
# pop_p


###
### (2) Build the model ----
###

# Subset dataframes to individuals only in the kinship matrix ----
df_pl_sub <- filter(df_pl, genotype %in% rownames(K))
df_pl_sub$genotype <- droplevels(df_pl_sub$genotype, exclude = "BL")

df_lf_sub <- filter(df_lf, genotype %in% rownames(K))
df_lf_sub$genotype <- droplevels(df_lf_sub$genotype)

###
### 1.  Single year Model ----
###
str(df_pl_sub)

trait = "relgrowth" # colname for trait
df = df_pl_sub

## Leaf_trait ~ Population + Genotype + Ploidy + Block
## The way our data is coded, genotype is implicitly nested in population
scols <- c(trait, "ploidy", "population", "block", "genotype")
trait_wide <- df %>%
  dplyr::select(all_of(scols))
trait_wide <- trait_wide[complete.cases(trait_wide),]
dim(trait_wide)

# mod.full <- relmatLmer(as.formula(paste0(trait, "~ ploidy + (1|population) + (1|ID) + (1|block)")),
#                        trait_wide,
#                        relmat = list(ID = k)) # relation matrix

mod.full <- mmer(as.formula(paste0(trait, " ~ ploidy")),
                 random = ~ population + block + 
                   vsr(genotype, Gu=K),
                 rcov= ~ units, 
                 nIters = 1000,
                 data = trait_wide,
                 verbose = FALSE)

mod.full$convergence


## > Check model assumptions ----
# Normality of residuals (less important
# Make a histogram of residuals
#hist(resid(mod.full))
hist(mod.full$residuals)

# Homogeneity of variances
plot(mod.full)# starry night plot

r <- residuals(mod.full)
plot(mod.full$residuals ~ r$ploidy)
plot(mod.full$residuals ~ r$block)
# plot(mod.full$residuals ~ trait_wide$population[as.numeric(names(mod.full$residuals))])
plot(mod.full$residuals ~ r$population)

# > Parameter estimates ----

## Fixed effects
# Significance and df

## Estimate & SE
summary(mod.full) # also give SE

# Means & Confidence intervals
means <- emmeans(mod.full, specs = 'ploidy')
summary(means)

eff <- contrast(means,'pairwise')
# effects_summary = summary(eff,level = 1-0.05/15,infer = c(T,T)) # bonferroni correction of CI
# effects_summary$p.value = pmin(1,effects_summary$p.value*15)

summary(eff, infer = T) # estimate same as summary estimate

## Random effects
# Test significance w/ LRT
### Population
mod.pop <- mmer(as.formula(paste0(trait, " ~ ploidy")),
                random = ~ block +
                  vsr(genotype, Gu=K),
                rcov= ~ units, 
                nIters = 1000,
                data = trait_wide,
                verbose = FALSE)
mod.pop$convergence
anova(mod.full, mod.pop)

### Block
mod.block <- mmer(as.formula(paste0(trait, " ~ ploidy")),
                  random = ~ population +
                    vsr(genotype, Gu=K),
                  rcov= ~ units, 
                  nIters = 1000,
                  data = trait_wide,
                  verbose = FALSE)
mod.block$convergence
anova(mod.full, mod.block)

### Genotype
mod.geno <- mmer(as.formula(paste0(trait, " ~ ploidy")),
                 random = ~ population + block,
                 rcov= ~ units, 
                 nIters = 1000,
                 data = trait_wide,
                 verbose = FALSE)
mod.geno$convergence
anova(mod.full, mod.geno)

# 2.  Multi-year traits ----
str(df_pl_sub)

trait = "height" # colname for trait
df = df_pl_sub

scols <- c(paste0(trait, "yr1"), paste0(trait, "yr2"), "ploidy", "population", "block", "genotype")
trait_wide <- df %>%
  dplyr::select(all_of(scols)) %>%
  tidyr::pivot_longer(cols = c(paste0(trait, "yr1"), paste0(trait, "yr2")),
               names_to = "year",
               values_to = trait)
trait_wide <- trait_wide[complete.cases(trait_wide),]
trait_wide$year <- as.factor(trait_wide$year)
dim(trait_wide)

### Leaf_trait ~ Population + Genotype + Ploidy + Block + Year
# mod.full <- relmatLmer(as.formula(paste0(trait, "~ ploidy + (1|population) + (1|ID) + (1|block) + (1|year)")),
#                        trait_wide,
#                        relmat = list(ID = k)) # relation matrix

mod.full <- mmer(as.formula(paste0(trait, " ~ ploidy")),
                 random = ~ population + block + year +
                   vsr(genotype, Gu=K),
                 rcov= ~ units, 
                 nIters = 10,
                 data = trait_wide,
                 verbose = FALSE)
mod.full$convergence

# mod.red <- mmer(as.formula(paste0(trait, " ~ 1")),
#                  random = ~ population + block + year +
#                    vsr(genotype, Gu=K),
#                  rcov= ~ units, 
#                  nIters = 10,
#                  data = trait_wide,
#                  verbose = FALSE)
# mod.red$convergence

# summary(mod.full)

# > Check model assumptions ----
# Normality of residuals (less important)
# Make a histogram of residuals
#hist(resid(mod.full))
hist(mod.full$residuals)

# Homogeneity of variances
plot(mod.full)# starry night plot

r <- residuals(mod.full)
plot(mod.full$residuals ~ r$ploidy)
plot(mod.full$residuals ~ r$block)
# plot(mod.full$residuals ~ trait_wide$population[as.numeric(names(mod.full$residuals))])
plot(mod.full$residuals ~ r$population)
plot(mod.full$residuals ~ r$year)

# > Parameter estimates ----
#  Fixed effects
## Significance and df
# car::Anova(mod.full, type = "3") # LRT
# anova(mod.full) # anova
# anova(mod.full, mod.red) # LRT

means_ploidy <- emmeans(mod.full, specs = 'ploidy')
eff_ploidy <- contrast(means_ploidy,'pairwise')
summary(eff_ploidy, infer = T) # estimate same as summary estimate

## Estimate & SE
summary(mod.full) # also give SE

# Random effects
## Significance w/ LRT
# lmerTest::ranova(mod.full)

### test population
mod.pop <- mmer(as.formula(paste0(trait, " ~ ploidy")),
                 random = ~ block + year +
                   vsr(genotype, Gu=K),
                 rcov= ~ units, 
                 nIters = 10,
                 data = trait_wide,
                 verbose = FALSE)
mod.pop$convergence
anova(mod.full, mod.pop)

### test block
mod.block <- mmer(as.formula(paste0(trait, " ~ ploidy")),
                 random = ~ population + year +
                   vsr(genotype, Gu=K),
                 rcov= ~ units, 
                 nIters = 20,
                 data = trait_wide,
                 verbose = FALSE)
mod.block$convergence
anova(mod.full, mod.block)

### test year
mod.year <- mmer(as.formula(paste0(trait, " ~ ploidy")),
                 random = ~ population + block +
                   vsr(genotype, Gu=K),
                 rcov= ~ units, 
                 nIters = 20,
                 data = trait_wide,
                 verbose = FALSE)
mod.year$convergence
anova(mod.full, mod.year)

### test genotype - drop everything with the kinship mat?
mod.geno <- mmer(as.formula(paste0(trait, " ~ ploidy")),
                 random = ~ population + block + year,
                 rcov= ~ units, 
                 nIters = 10,
                 data = trait_wide,
                 verbose = FALSE)
mod.geno$convergence
anova(mod.full, mod.geno)
# an.g <- anova(mod.full, mod.geno)
# an.g$Chisq[2]
# gsub(pattern = "[^0-9.e-]", replacement = "", x = an.g$PrChisq[2])

# Variances 
# sum <- as.data.frame(VarCorr(mod.full)) # variances for random effects
# sum[,c(1,4)] # grab sd
# summary(mod.full)

# R2 // check fit ??? XX
# marginal R squared - variance explained by fixed effects only
# conditional R squared - variance explained by fixed and random effects
# MuMIn::r.squaredGLMM(mod.full)

# > Heritability ----
var <- summary(mod.full)$varcomp
h2 <- vpredict(mod.full, h2 ~ V4 / (V3 + V4 + V1 + V2 + V5)) # narrow sense = geno / (all + error)
h2


# (3) Environment ----

# Load climate data ----

##  Downloaded from ClimateNA ----
clm <- read.csv("~/cg_andropogon/scripts/Climate Data - ClimateNA/androclimatedata_popcodes.csv")

# Downloaded growing degree days > 10 degrees C from DAYMET
GDD <- read.csv("~/cg_andropogon/scripts/GDD Calculations/andro_yearly_gdd.csv")
clm$GDD <- colMeans(GDD[,2:17]) # average from 1980 to 2000

bou <- colMeans(clm[1:2,3:31]) # take the average environment from each boulder population
clm_df <- clm[2:16,]
clm_df[1,3:31] <- bou 

# PCA on env variables ----
dim(clm_df)
pc <- prcomp(clm_df[,3:31], center = T, scale = F)

# PC1 + cor: MAR, MCMT, MAT, EMT, Eref, FFP, eFFP, NFFD, DD1040, DD5, GDD, DD18, EXT, MWMT
# PC1 - cor: DD_18, lat, DD_0, TD, PAS, bFFP
# PC15: flipped PC1
# PC3 + cor: MSP, MAP, long, RH, CMI
# PC3 - cor: CMD, AHM, SHM, MAR

# pc_var <- summary(pc)
# pc_var$importance[2,]*100

# PC1 = 91.665%, PC3 = 2.721%, PC15 = 0.00%

# Correlation between env var
# clm_cor <- cor(clm_df[,c(3:4, 6:31)]) # exclude elevation (el)
# pheatmap::pheatmap(clm_cor)


# MAP is strongly positively correlated with CMI, MSP, long, and RH
# MAT and GDD are strongly positively correlated with MCMT, Eref (evapotranspiration), MWMT, NFFD, EMT, DD18, DD1040, DD5, eFFP, EXT, and FFP 
# AHM is strongly positively correlated with MAR, CMD, and SHM (negatively correlated with MAP group)
# TD (continentality) is strongly positively correlated with DD_0, PAS, lat, bFFP, DD_18 (parameters dealing with freezing) (negatively corrlated with MAT group)

# PC Loadings ----
pheatmap::pheatmap(pc$rotation[,1:3]*100)
round(pc$rotation[,1]*100,3) %>% sort()  # PC1
round(pc$rotation[,3]*100,3) %>% sort() # PC3

ggsave

# PC1 = + GDD, DD5, DD1040; - DD_0, DD_18
# PC2 = + GDD>10, DD>5, DD>10<40, DD<0, DD<18; - CMD
# PC3 = - CMD, + MAP (CMD =  (mm)

# Select env & PCs var for hyp testing ----
## Mean annual temp (MAT)
## Mean annual precipitation (MAP)
## Growing degree days over 10C (GDD)
## Annual heat-moisture index (AHM)
## Continentality (TD)
## PC1 & PC3 captures all the above
clim_var <- dplyr::select(clm_df, ID2, MAT, MAP, GDD, AHM, TD) # GDD and MAT are redundant
colnames(clim_var) <- c("pop_short", "MAT", "MAP", "GDD", "AHM", "TD")

clim_var$PC1 <- pc$x[,1]
clim_var$PC3 <- pc$x[,3]
# pheatmap::pheatmap(cor(clim_var[,2:8]))
pheatmap::pheatmap(cor(clm_df[,c(3:4,6:31)], pc$x[,c(1,3)])) # supplemental fig

head(clim_var)
str(clim_var)

# Export final climate variable set
# write.csv(clim_var, file = paste0("~/cg_andropogon/data/final_env_var.", Sys.Date(), ".csv"), quote = F)

#  Subset dataframes to ind only in the kinship matrix and add clim vars to df ----
df_pl_tmp <- filter(df_pl, genotype %in% rownames(K))
df_pl_tmp$genotype <- droplevels(df_pl_tmp$genotype, exclude = "BL")
df_pl_clim <- merge(df_pl_tmp, clim_var, "pop_short", all.x = T, sort = T)
str(df_pl_clim)

df_lf_tmp <- filter(df_lf, genotype %in% rownames(K))
df_lf_tmp$genotype <- droplevels(df_lf_sub$genotype)
df_lf_clim <- merge(df_lf_tmp, clim_var, "pop_short", all.x = T, sort = T)
str(df_lf_clim)

###
### 1. Single year traits ----
###

df_lf_clim$genotype <- droplevels(df_lf_clim$genotype)

trait = "stomdenyr2" # colname for trait
df = df_lf_clim

str(df_lf_clim)

scols <- c(trait, "ploidy", "population", "block", "genotype", "PC1", "PC3")
# scols <- c(trait, "ploidy", "population", "block", "genotype", "MAT")
trait_wide <- df %>%
  dplyr::select(all_of(scols))
trait_wide <- trait_wide[complete.cases(trait_wide),]
head(trait_wide)
# str(trait_wide)

# > Build models ----
mod.clim <- mmer(as.formula(paste0(trait, " ~ ploidy + PC1 + PC3")),
                 random = ~ population + block +
                   vsr(genotype, Gu=K),
                 rcov= ~ units,
                 nIters = 10,
                 data = trait_wide,
                 verbose = FALSE)

mod.clim$convergence
summary(mod.clim)

# > Model assumptions ----
hist(mod.clim$residuals)
plot(mod.clim)# starry night plot

r <- residuals(mod.clim)
plot(mod.clim$residuals ~ r$ploidy)
plot(mod.clim$residuals ~ r$block)
plot(mod.clim$residuals ~ r$population)
plot(mod.clim$residuals ~ r$PC1)
plot(mod.clim$residuals ~ r$PC3)

# > Test fixed effects ----

anova(mod.clim)

means_ploidy <- emmeans(mod.clim, specs = 'ploidy')
eff_ploidy <- contrast(means_ploidy,'pairwise')
summary(eff_ploidy, infer = T) # estimate same as summary estimate
  
# > Test random effects ----
# Population no longer represents site but instead the non-genetic and non-climate effects of population
mod.pop.clim <- mmer(as.formula(paste0(trait, " ~ ploidy + PC1 + PC3")),
                random = ~ block +
                  vsr(genotype, Gu=K),
                rcov= ~ units, 
                nIters = 10,
                data = trait_wide,
                verbose = FALSE)
mod.pop.clim$convergence
anova(mod.clim, mod.pop.clim)

mod.block.clim <- mmer(as.formula(paste0(trait, " ~ ploidy + PC1 + PC3")),
                  random = ~ population +
                    vsr(genotype, Gu=K),
                  rcov= ~ units, 
                  nIters = 10,
                  data = trait_wide,
                  verbose = FALSE)
mod.block.clim$convergence
anova(mod.clim, mod.block.clim)

mod.geno.clim <- mmer(as.formula(paste0(trait, " ~ ploidy + PC1 + PC3")),
                 random = ~ population + block,
                 rcov= ~ units, 
                 nIters = 10,
                 data = trait_wide,
                 verbose = FALSE)
mod.geno.clim$convergence
anova(mod.clim, mod.geno.clim)


###
### 2. Multi-year traits ----
###
str(df_pl_clim)

trait = "preretillersEXP" # colname for trait
df = df_pl_clim

scols <- c(paste0(trait, "yr1"), paste0(trait, "yr2"), "ploidy", "population", "block", "genotype", "PC1", "PC3")
trait_wide <- df %>%
  dplyr::select(all_of(scols)) %>%
  tidyr::pivot_longer(cols = c(paste0(trait, "yr1"), paste0(trait, "yr2")),
                      names_to = "year",
                      values_to = trait)
trait_wide <- trait_wide[complete.cases(trait_wide),]
trait_wide$year <- as.factor(trait_wide$year)
head(trait_wide)

# > Build model ----
mod.clim <- mmer(as.formula(paste0(trait, " ~ ploidy + PC1 + PC3")),
                 random = ~ population + block + year +
                   vsr(genotype, Gu=K),
                 rcov= ~ units, 
                 nIters = 10,
                 data = trait_wide,
                 verbose = FALSE)
mod.clim$convergence
summary(mod.clim)

# > Model assumptions ----
hist(mod.clim$residuals)
plot(mod.clim)# starry night plot

r <- residuals(mod.clim)
plot(mod.clim$residuals ~ r$ploidy)
plot(mod.clim$residuals ~ r$block)
plot(mod.clim$residuals ~ r$population)
plot(mod.clim$residuals ~ r$PC1)
plot(mod.clim$residuals ~ r$PC3)

# > Test fixed effects ----
anova(mod.clim)

means_ploidy <- emmeans(mod.clim, specs = 'ploidy')
eff_ploidy <- contrast(means_ploidy,'pairwise')
summary(eff_ploidy, infer = T) # estimate same as summary estimate

# > Test random effects ----
# Population no longer represents site but instead the non-genetic and non-climate effects of population
mod.pop.clim <- mmer(as.formula(paste0(trait, " ~ ploidy + PC1 + PC3")),
                     random = ~ block + year +
                       vsr(genotype, Gu=K),
                     rcov= ~ units, 
                     nIters = 10,
                     data = trait_wide,
                     verbose = FALSE)
mod.pop.clim$convergence
anova(mod.clim, mod.pop.clim)

mod.block.clim <- mmer(as.formula(paste0(trait, " ~ ploidy + PC1 + PC3")),
                       random = ~ population + year +
                         vsr(genotype, Gu=K),
                       rcov= ~ units, 
                       nIters = 10,
                       data = trait_wide,
                       verbose = FALSE)
mod.block.clim$convergence
anova(mod.clim, mod.block.clim)

mod.year.clim <- mmer(as.formula(paste0(trait, " ~ ploidy + PC1 + PC3")),
                       random = ~ population + block +
                         vsr(genotype, Gu=K),
                       rcov= ~ units, 
                       nIters = 10,
                       data = trait_wide,
                       verbose = FALSE)
mod.year.clim$convergence
anova(mod.clim, mod.year.clim)

mod.geno.clim <- mmer(as.formula(paste0(trait, " ~ ploidy + PC1 + PC3")),
                      random = ~ population + block + year,
                      rcov= ~ units, 
                      nIters = 10,
                      data = trait_wide,
                      verbose = FALSE)
mod.geno.clim$convergence
anova(mod.clim, mod.geno.clim)

# (4) Plot traits vs Env ----

plot_env <- function(df, y_var){
  genotype_means <- aggregate(as.formula(paste0(y_var, "~ pop_short + genotype + ploidy + PC1 + PC3")), df, FUN = mean)
  
  p1 <- ggplot(genotype_means,  
               aes(x = reorder(.data[["pop_short"]], .data[["PC1"]]), 
                   y = .data[[y_var]], 
                   fill = ploidy)) + 
    geom_boxplot()+
    geom_point(position = position_jitterdodge(jitter.width = 0)) +
    xlab("Population ordered by ENV PC1") +
    ylab(paste0(y_var)) +
    theme_bw() +
    scale_fill_manual(values = c("white", "darkgray"))
  
  p3 <- ggplot(genotype_means,  
               aes(x = reorder(.data[["pop_short"]], .data[["PC3"]]), 
                   y = .data[[y_var]], 
                   fill = ploidy)) + 
    geom_boxplot()+
    geom_point(position = position_jitterdodge(jitter.width = 0)) +
    xlab("Population ordered by ENV PC3") +
    ylab(paste0(y_var)) +
    theme_bw() +
    scale_fill_manual(values = c("white", "darkgray"))

ggarrange(p1, p3, 
          ncol = 1,
          nrow = 2,
          heights = c(2.25, 2.25)
)
  
ggsave(paste0("~/cg_andropogon/data_explore/",y_var, ".environment.",Sys.Date(),".pdf"),
       height = 6,
       width = 8)
}

## 1. All traits vs env ----
## plots for whole-plant traits
# for (i in colnames(df_pl)[8:24]){
#   plot_env(df_pl_clim, i)
# }

## plots for leaf traits
# for (i in colnames(df_lf)[8:21]){
#   plot_env(df_lf_clim, i)
# }

# plot_env(df_lf_clim, "stomlenyr2" )

## 2. Figures for publication ----
genotype_means <- aggregate(stomdenyr2 ~ pop_short + genotype + ploidy + PC1 + PC3, df_lf_clim, FUN = mean)
hex_colors <- c("#FF0000", "#00FF00", "#0066FF", "#FFFF00", "#00FFFF", "#FF00FF", "#FFA500", "#800080", "#008000", "#FF6347", "#00CED1", "#9A2FE2", "#7FFFD4", "#FFD700")
# y_name <- expression(Stomatal ~ density ~ (m^-2))

## * Stomata vs PC3 ----

# Estimate uncertainty as sqrt(var(intercept) + PC3^2*var(slope))
mod.clim$Beta
mod.clim$VarBeta
intercept <- mod.clim$Beta[1,3]
var.intercept <- mod.clim$VarBeta[1,1]
var.PC3 <- mod.clim$VarBeta[4,4] # Variance of PC3
slope.PC3 <- mod.clim$Beta[4,3] 

pc3.values <- seq(-645, 385, 10)
se.PC3 <- sqrt(var.intercept + (pc3.values^2)*var.PC3) # SE at each point on the x-axis

predicted_values <- slope.PC3 * pc3.values + intercept 

upper.se <- 2 * se.PC3 + predicted_values
lower.se <- predicted_values - 2* se.PC3

data.frame(upper.se, pc3.values)

# Plot

pc3_p <- ggplot(genotype_means,  aes(x = PC3, y = stomdenyr2)) +
  geom_point(aes(color = pop_short, shape = ploidy), 
             size = 3, alpha = 0.6) +
  # geom_smooth(method = "lm")
  geom_abline(mapping = aes(intercept = intercept, slope = slope.PC3), # values from clim model
              color = "black", linewidth = 1) +
  geom_line(data=data.frame(upper.se, pc3.values),
  aes(x = pc3.values, y = upper.se), linetype = "dashed") +
  geom_line(data=data.frame(lower.se, pc3.values),
  aes(x = pc3.values, y = lower.se), linetype = "dashed") +
  # geom_abline(mapping = aes(intercept = 241.1392 + 2* 7.880465, slope = slope.PC3), # values from clim model
  #             color = "black", size = .7, linetype="dashed") +
  # geom_abline(mapping = aes(intercept = 241.1392 - 2 * 7.880465, slope = slope.PC3), # values from clim model
  #             color = "black", size = .7, linetype="dashed") +
  xlab("PC3") +
  ylab(expression('Stomata density (m'^2*')')) +
  theme_bw() +
  scale_color_manual(values = hex_colors, name = "Population") +
  guides(shape = 'none')+
  theme(text = element_text(size = 12))

pc3_p

ggsave(pc3_p, file = paste0("~/cg_andropogon/Figures/PC3vsStomataDensity_twoSEs.", Sys.Date(), ".jpeg"),
       width = 6,
       height = 5,
       unit = "in")

## * Plot PCs vs climate variables ----
clim_var$pop_short <- as.factor(clim_var$pop_short)
clim_p <- clim_var %>%
  filter(!pop_short %in% c("Site")) %>%
  ggplot(aes(x = PC3, y = MAP)) +
  geom_smooth(method = "lm", col = "black") +
  geom_point(aes(color = pop_short), size = 3) +
  theme_bw() +
  scale_color_manual(values = hex_colors, name = "Population") +
  theme(text = element_text(size = 12)) +
    ylab("MAP (mm)") +
  guides(color = "none")

clim_p

ggsave(clim_p, file = paste0("~/cg_andropogon/Figures/PC3vsMAP.", Sys.Date(), ".jpeg"),
       width = 5,
       height = 5,
       unit = "in")

# ggarrange(pc1_p, pc3_p, nrow = 1, ncol = 2)

ggarrange(clim_p, pc3_p,stom_p, nrow = 2, ncol = 2)

ggsave(paste0("~/cg_andropogon/data_explore/PC3vsMAP.", Sys.Date(), ".jpeg"),
       width = 6,
       height = 4,
       unit = "in")


## test test
rg <- qdrg(~ ploidy + PC1 + PC3, data = mod.clim$data, coef = mod.clim$Beta$Estimate, 
           vcov = mod.clim$VarBeta) # create reference grid




# (5) Consistency of stomata impressions ----
## 1. Load second and third impression data ----
## Already in long format
impressions <- read.csv("~/cg_andropogon/data/2022/2022 stomata data - Stomata 2-3_d03142024.csv")
head(impressions)
str(impressions)

impressions$stomlenyr2 <- rowMeans(impressions[,5:9], na.rm = T) # calculate mean length
impressions$stomdenyr2 <- (impressions$stom_num/(impressions$area_um2/1000000)) # density = # stomata / area (mm^2)
impressions$stomspiyr2 <- impressions$stomdenyr2/sqrt(impressions$stomlenyr2) # SPI = d/sqrt(mean length)

impressions$geno_leaf_imp <- paste0(impressions$genotype, ".", impressions$leaf, ".",impressions$impression)
impressions$geno_leaf <- paste0(impressions$genotype, ".", impressions$leaf)
str(impressions)

## 2. Merge dataframes ----
stom_df <- dplyr::select(df_lf_sub, plant, leaf, stomlenyr2, stomdenyr2, stomspiyr2)
stom_df$impression <- 1
# stom_df$geno_leaf_imp <- paste0(stom_df$plant, ".", stom_df$leaf, ".",stom_df$impression) %>% as.factor()
# stom_df$geno_leaf <- paste0(stom_df$plant, ".", stom_df$leaf)
str(stom_df)

imp_sub <- impressions[,c(1:3,12:14)]
imp_sub <- imp_sub[complete.cases(imp_sub),]
imp_sub <- imp_sub %>% relocate(impression, .after = last_col())

colnames(stom_df) <- colnames(imp_sub)
colnames(stom_df)
colnames(imp_sub)

# Subset to genotypes with second and third impressions
stom_df_sub <- stom_df[stom_df$genotype %in% imp_sub$genotype,]

# Stack two datasets
all_stom <- rbind(stom_df_sub, imp_sub)
# all_stom <- all_stom[complete.cases(all_stom),]
dim(all_stom)
str(all_stom)



# Add metadata 
meta_sub <- df_lf_sub[df_lf_sub$plant %in% all_stom$genotype, c(1:2, 4:6)] %>% unique()
# colnames(meta_sub)[1] <- "genotype"
colnames(all_stom)
colnames(meta_sub)
all_stom_df <- merge(all_stom, meta_sub, by.x = "genotype", by.y = "plant")
colnames(all_stom_df)[c(1,7)] <- c("plant", "genotype")

all_stom_df$impression <- as.factor(all_stom_df$impression)
all_stom_df$genotype <- as.factor(all_stom_df$genotype)
all_stom_df$genotype <- droplevels(all_stom_df$genotype)

str(all_stom_df)

# drop outliers
all_stom_df <- all_stom_df[all_stom_df$stomdenyr2 < 1000,]

# 3. Plot ----
# Take means of leaves
impr_means <- aggregate(stomspiyr2 ~ impression + plant + genotype + ploidy + pop_short, all_stom_df, mean )

# Overall distributions are similar across impressions
# ggplot(impr_means, aes(x = impression, y = stomspiyr2)) +
#   geom_boxplot() +
#   geom_jitter(alpha = 0.1)
# 
# ggplot(impr_means, aes(x = impression, y = stomspiyr2)) +
#   geom_violin()

# By population
# ggplot(impr_means, aes(x = impression, y = stomspiyr2)) +
#   geom_boxplot() +
#   geom_jitter(alpha = 0.1) +
#   facet_wrap(~ pop_short)

# By ploidy (supplemental fig)
axis_names <- c( expression(paste0("6", italic("x"))),
                 expression(paste0("9", italic("x"))) )

spi_p <- ggplot(impr_means, aes(x = impression, y = stomspiyr2, fill = ploidy)) +
  geom_boxplot() +
  geom_point(aes(fill = ploidy), alpha = 0.2, position = position_jitterdodge(jitter.width = 0.1)) +
  # geom_jitter(alpha = 0.1, width = 0.1) +
  # facet_wrap(~ ploidy) +
  theme_bw() +
  xlab("Impression") + ylab("Mean SPI") +
  scale_fill_manual(values = c("white", "#666666")) +
  guides(fill = "none")

den_p <- ggplot(impr_means, aes(x = impression, y = stomdenyr2, fill = ploidy)) +
  geom_boxplot() +
  geom_point(aes(fill = ploidy), alpha = 0.2, position = position_jitterdodge(jitter.width = 0.1)) +
  # geom_jitter(alpha = 0.1, width = 0.1) +
  # facet_wrap(~ ploidy) +
  theme_bw() +
  xlab("Impression") + ylab(expression('Mean stomata density (m'^2*')')) +
  scale_fill_manual(values = c("white", "#666666")) +
  guides(fill = "none")

len_p <- ggplot(impr_means, aes(x = impression, y = stomlenyr2, fill = ploidy)) +
  geom_boxplot() +
  # geom_jitter(alpha = 0.1, width = 0.1) +
  geom_point(aes(fill = ploidy), alpha = 0.2, position = position_jitterdodge(jitter.width = 0.1)) + 
  # facet_wrap(~ ploidy) +
  theme_bw() +
  xlab("Impression") + ylab(expression(paste("Mean stomata length (", mu, "m)"))) +
  scale_fill_manual(values = c("white", "#666666")) +
  guides(fill = "none")

impr_p <- ggarrange(spi_p, den_p, len_p, nrow = 1, ncol = 3)
ggsave(impr_p, file = paste0("~/cg_andropogon/Figures/stomataimpressions_by_ploidycolor.", Sys.Date(),".jpeg"),
       height = 4, width = 12, unit = "in")

# By plant
ggplot(impr_means , aes(x = impression, y = stomspiyr2, 
                        color = pop_short,
                        shape = ploidy)) +
  geom_point() +
  geom_line(aes(group = plant, 
                color = pop_short
  ))+
  xlab('Impression') + 
  ylab('Mean SPI') +
  theme_minimal() +
  guides(color = guide_legend(title = "Population"),
         shape = guide_legend(title = "Ploidy"))

# Only impression 1 and 3
filter(impr_means, impression %in% c('2', '3')) %>% 
ggplot(aes(x = impression, y = stomspiyr2, 
                        color = pop_short,
                        shape = ploidy)) +
  geom_point() +
  geom_line(aes(group = plant, 
                color = pop_short
  ))+
  xlab('Impression') + 
  ylab('Mean stomata length') +
  theme_minimal() +
  guides(color = guide_legend(title = "Population"),
         shape = guide_legend(title = "Ploidy"))

# 4. Model effect of impression ----
# Can't use a linear mixed model to test the difference betwwen impressions as 

K_stom <- k[,colnames(k) %in% unique(all_stom_df$genotype)]
K_stom <- K_stom[rownames(K_stom) %in% unique(all_stom_df$genotype),]
dim(K_stom)
K_stom[1:5,1:5]
K_stom <- data.matrix(K_stom)

table(all_stom$genotype, all_stom$impression, exclude = T)
table(all_stom_df$pop_short, all_stom_df$impression, exclude = T)
table(all_stom_df$leaf, all_stom_df$impression, exclude = T)

all_stom_df$plant_leaf <- paste0(all_stom_df$plant,".",all_stom_df$leaf) %>% as.factor() # make uniquely nested leaves

mod.stom <- mmer(stomspiyr2 ~ ploidy,
                 random = ~ pop_short + block + plant_leaf + impression +
                 vsr(genotype, Gu=K_stom),
                 rcov= ~ units, 
                 nIters = 1000,
                 data = all_stom_df,
                 verbose = FALSE)

mod.stom$convergence
summary(mod.stom)

# Test importance of impression
mod.stom.impr <- mmer(stomspiyr2 ~ ploidy,
                      random = ~ pop_short + block + plant_leaf  +
                        vsr(genotype, Gu=K_stom),
                     rcov= ~ units, 
                     nIters = 1000,
                     data = all_stom_df,
                     verbose = FALSE)
mod.stom.impr$convergence
anova(mod.stom, mod.stom.impr) # Impression contributes to a significant proportion of variation

## Calculate means for just the three impressions (for figure)
mod.impr <- lm(stomlenyr2 ~ impression, all_stom_df)
impr.means <- emmeans(mod.impr, specs = 'impression')
impr.diffs <- contrast(impr.means, 'pairwise')
impr.diffs

'The third impression is signficantly different than the first two impressions in SPI (3 SPI unit) and density (14 and 22 stom/area). The first impression is signficantly different from the last two in length(2 um) .'

# Is the effect of ploidy consistent across impressions? Subset to third impression
third_stom_df <- all_stom_df %>% filter(impression == '3')
K_3stom <- k[,colnames(k) %in% unique(third_stom_df$genotype)]
K_3stom <- K_3stom[rownames(K_3stom) %in% unique(third_stom_df$genotype),]
dim(K_3stom)
K_3stom[1:5,1:5]
K_3stom <- data.matrix(K_3stom)

mod.3stom <- mmer(stomspiyr2 ~ ploidy,
                 random = ~ pop_short + block + 
                   vsr(genotype, Gu = K_3stom),
                 rcov= ~ units, 
                 nIters = 1000,
                 data = third_stom_df,
                 verbose = FALSE)

mod.3stom$convergence

anova(mod.3stom)

stom3.means <- emmeans(mod.3stom, specs = 'ploidy')
summary(stom3.means)
stom3.eff <- contrast(stom3.means,'pairwise')
stom3.eff
