### Title: Common garden MANOVA and RDA
### Author: Alyssa Phillips
### Date: 4/3/23

library(ggplot2)
library(dplyr)
# install_github("vqv/ggbiplot")
library(ggbiplot)
library(ggpubr)

# Adapted from scripts by Grey Monroe (UC Davis)

# (1) Load data ----

## Trait data ----
df_pl <- read.csv("~/cg_andropogon/data/plant_traits.2024-03-27.csv", header = T)
df_pl[,1:9] <- lapply(df_pl[,1:9], as.factor)
str(df_pl)

df_lf <- read.csv("~/cg_andropogon/data/leaf_traits.2024-03-27.csv", header = T)
df_lf[,1:7] <- lapply(df_lf[,1:7], as.factor)
str(df_lf)

## Kinship matrix ----
k <- read.table("~/Andropogon/kinship/all.andro.lowcov.100k.miss20.centereddiag.2023-09-27.txt")
colnames(k) <- rownames(k)
k <- as.matrix(k)
k[1:5,1:5]
dim(k)

# Subset k
K <- k[,colnames(k) %in% unique(df_lf$genotype)]
K <- K[rownames(K) %in% unique(df_lf$genotype),]
dim(K)
K[1:5,1:5]

# Subset dataframes
df_pl_sub <- filter(df_pl, genotype %in% rownames(K))
df_pl_sub$genotype <- droplevels(df_pl_sub$genotype, exclude = "BL")

df_lf_sub <- filter(df_lf, genotype %in% rownames(K))
df_lf_sub$genotype <- droplevels(df_lf_sub$genotype)

# (2) Estimate BLUPs for each trait ----
# df <- dplyr::select(pl, plant, genotype, pop_short, block, heightyr1)
# df <- df[complete.cases(df),]
# df[,1:4] <- lapply(df[,1:4], factor)
# df

## Make output matrix ----
# Genotypes values will be output alphabetically 
# rownames(blups) == sort(rownames(K)) # TRUE
blup_mat <- matrix(NA, nrow = 78, ncol = 16) # 78 genos by (15 traits + 1 geno) cols
blup_df <- as.data.frame(blup_mat)
colnames(blup_df) <- c("genotype", "basalgrowthLOG", "relgrowth", "tillgrowth", "heightgrowth",
                       "stomdenyr2", "stomlenyr2", "stomspiyr2", "ldmcyr2",
                       "height", "tillers", "preretillersEXP",
                       "slaINT", "length", "width", "thick")

blup_df[,"genotype"] <- sort(rownames(K))

## Single year traits ----
# single_year <- c("basalgrowthLOG", "relgrowth", "tillgrowth", "heightgrowth") # df_pl
single_year <- c("stomdenyr2", "stomlenyr2", "stomspiyr2", "ldmcyr2") # df_lf

df = df_lf_sub
# t = as.character("relgrowth")

for (t in single_year){
  mod.blup <- mmer(as.formula(paste0(t, "~ 1 + ploidy")),
                   random = ~vsr(genotype, Gu = K) + pop_short + block ,
                   rcov= ~units, 
                   nIters = 1000,
                   data = df, verbose = FALSE)
  
  print(mod.blup$convergence) # check convergence
  
  blup_df[,t] <- as.data.frame(mod.blup$U$`u:genotype`)
  print(t)
}

## Multi-year traits ----
multi_year <- c("height", "tillers", "preretillersEXP") # df_pl
# multi_year <- c("slaINT", "length", "width", "thick") # df_lf

df = df_pl_sub
# t = as.character("slaINT")
for (t in multi_year){
  scols <- c(paste0(t, "yr1"), paste0(t, "yr2"), "ploidy", "pop_short", "block", "genotype")
  trait_wide <- df %>%
    dplyr::select(all_of(scols)) %>%
    tidyr::pivot_longer(cols = c(paste0(t, "yr1"), paste0(t, "yr2")),
                        names_to = "year",
                        values_to = t)
  trait_wide <- trait_wide[complete.cases(trait_wide),]
  trait_wide$year <- as.factor(trait_wide$year)
  
  mod.blup <- mmer(as.formula(paste0(t, "~ 1 + ploidy")),
                   random = ~vsr(genotype, Gu = K) + pop_short + block + year,
                   rcov= ~units, 
                   nIters = 1000,
                   data = trait_wide, verbose = FALSE)
  
  print(mod.blup$convergence) # check convergence
  
  blup_df[,t] <- as.data.frame(mod.blup$U$`u:genotype`)
  print(t)
}

# Export blups ----
# write.csv(blup_df,paste0("~/cg_andropogon/data/blups.", Sys.Date(), ".csv"), row.names = F, quote = F)

# Reload blups ----
blup_df <- read.csv("~/cg_andropogon/data/blups.2024-03-27.csv")
str(blup_df)

## Look at correlation between estimated BLUPs and mean trait values
### Should be pretty close but not identical
# g <- blup_df$genotype %>% as.factor()
# df_means <- aggregate(tillgrowth ~ genotype + pop_short, FUN = mean, df_pl_sub) %>%
#   arrange(genotype)
# cor(blup_df[,"tillgrowth"], df_means[g, "tillgrowth"], use="complete")
# plot(blup_df[,"tillgrowth"], df_means[g, "tillgrowth"])


# (3) Calculate Z-scores ----
traits_scaled <- scale(blup_df[,c(2:5,8:16)]) # only one stomata trait
traits_scaled <- scale(blup_df[,c(2:4,7:16)])

## correlation among traits
# cor(blup_df[,2:16]) %>% pheatmap::pheatmap()

# (4) Run PCA ----
pca <- prcomp(traits_scaled)

# Make useable dataframe
g_meta <- df_pl_sub[df_pl_sub$genotype %in% blup_df$genotype, c("genotype","pop_short", "ploidy")] %>% 
  unique() %>%
  arrange(genotype)

## Load pop assignments
pop_assignment <- read.csv("~/Andropogon/structure/no_clones/pop_assignment2024-03-27.csv", header = T)

g_meta <- merge(g_meta, pop_assignment, by.x = "genotype", by.y = "geno", all.y = F, all.x = T)
g_meta$group <- as.factor(g_meta$group)
g_meta$group[is.na(g_meta$group)] <- "2" 
## 2 = West 
## 1 = East

# blup_df$genotype == g_meta$genotype # TRUE

pca_df <- data.frame(pca$x, g_meta)
head(pca_df)


## Run on just Year 2 trait means ----
# yr2_lf_traits <- c("widthyr2", "lengthyr2", "slayr2", "ldmcyr2", "thickyr2", "stomspiyr2")
# df_lf_yr2 <- dplyr::select(df_lf_sub, genotype, widthyr2, lengthyr2, slayr2, ldmcyr2, thickyr2, stomspiyr2)
# 
# yr2_pl_traits <- c("preretillersEXPyr2", "basalyr2", "tillersyr2", "heightyr2", "relgrowth")
# df_pl_yr2 <- dplyr::select(df_pl_sub, genotype, all_of(yr2_pl_traits))
# 
# num_traits <- length(yr2_lf_traits) + length(yr2_pl_traits)
# num_genos <- unique(df_lf_yr2$genotype) %>% length()
# trait_means <- matrix(NA, nrow = num_genos, ncol = num_traits +1) %>% data.frame()
# colnames(trait_means) <- c("genotype", yr2_lf_traits, yr2_pl_traits)
# trait_means$genotype <- aggregate(as.formula(paste0("basalyr2", "~ genotype")), df_pl_yr2, FUN = mean)[,1]
# 
# for (trait in yr2_lf_traits){
#   trait_means[, trait] <- aggregate(as.formula(paste0(trait, "~ genotype")), df_lf_yr2, FUN = mean)[,2]
# }
# 
# for (trait in yr2_pl_traits){
#   trait_means[, trait] <- aggregate(as.formula(paste0(trait, "~ genotype")), df_pl_yr2, FUN = mean)[,2]
# }
# 
# traits_scaled <- scale(trait_means[,2:12])
# 
# pca <- prcomp(traits_scaled)

# (5) Plot ----

pc <- summary(pca)$importance[2,]
hex_colors <- c("#FF0000", "#00FF00", "#0066FF", "#FFFF00", "#00FFFF", "#FF00FF", "#FFA500", "#800080", "#008000", "#FF6347", "#00CED1", "#9A2FE2", "#7FFFD4", "#FFD700")

ggbiplot(pca,
         choices = 1:2,
         groups = g_meta$group,
         ellipse = TRUE,
         ellipse.prob = 0.9,
         ellipse.fill = F,
         ellipse.linewidth = 0.75,
         var.axes = F,
         alpha = 0,
         scale = 0, # disable scaling
) +
  # scale_linetype_manual(values = c("dashed", "solid")) +
  geom_point(data = pca_df, aes(x = PC1, y = PC2, colour = pop_short, 
                                # shape = ploidy
                                )) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'vertical', 
        legend.position = 'right')  +
  theme_bw() +
  xlab(paste0("PC1 (", round(pc[1] * 100, 1), "%)")) +
  ylab(paste0("PC2 (", round(pc[2] * 100, 1), "%)")) +
  theme(
    legend.position = "none",
        panel.background = element_blank(),
        panel.grid = element_blank()) +
  scale_color_manual(values = c("black", "gray", hex_colors ),
                     # breaks = unique(pca_df$pop_short),
                     name = "Population") + # black == east
  # geom_line(aes(linetype = g_meta$group))
  guides(shape = "none", 
         color = guide_legend(override.aes = list(linetype = 0))) # remove line through dot

# PCA plot
# ggplot(pca_df, aes(x = PC1, y = PC2,
#                    col = pop_short,
#                    shape = ploidy,
#                    # col = group
#                    )) +
#   geom_point(size = 2) +
#   # scale_shape_manual(values=c(21,24)) +
#   # scale_color_manual(values = c("gray", "black")) +
#   scale_color_manual(values = hex_colors) +
#   theme_bw(base_size = 12) +
#   xlab(paste0("PC1 (", round(pc[1] * 100, 1), "%)")) +
#   ylab(paste0("PC2 (", round(pc[2] * 100, 1), "%)")) +
#   guides(shape = "none",
#          # fill = guide_legend("Population", override.aes = list(shape = 21)),
#          # col = "none",
#          col = guide_legend("Population")
#          )

ggsave(paste0("~/cg_andropogon/Figures/traits_pca.byPop.", Sys.Date(), ".jpeg"),
       units = "in",
       width = 43,
       height = 6)

# Biplot
trait_names <- c("Basal growth","Relative aboveground growth","Tiller growth","Height growth"   , "SPI","LDMC","Height","Tiller count","Percent flowering tillers", "SLA","Leaf length","Leaf width","Leaf thickness")
rownames(pca$rotation) <- trait_names

p_group <- ggbiplot(pca,
         choices = 1:2,
         groups = g_meta$group,
         ellipse = TRUE,
         ellipse.fill = F,
         ellipse.linewidth = 0.75,
         var.axes = F,
         # circle = TRUE,
         ellipse.prob = 0.9,
         # scale = 0,
         varname.size = 3
         ) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'vertical', 
        legend.position = 'right')  +
   theme_bw() +
  xlab(paste0("PC1 (", round(pc[1] * 100, 1), "%)")) +
  ylab(paste0("PC2 (", round(pc[2] * 100, 1), "%)")) +
  theme(legend.position = "none", 
        panel.background = element_blank(),
        panel.grid = element_blank()) +
  scale_color_manual(values = c("#E1BE6A","#40B0A6" )) # 2 is green

p_ploidy <- ggbiplot(pca,
                    choices = 1:2,
                    groups = g_meta$ploidy,
                    ellipse = TRUE,
                    # circle = TRUE,
                    ellipse.fill = F,
                    ellipse.linewidth = 0.75,
                    ellipse.prob = 0.9,
                    scale = 0,
                    var.axes = F,
                    # var.scale = 1
) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'vertical', 
        legend.position = 'right')  +
  theme_bw() +
  xlab(paste0("PC1 (", round(pc[1] * 100, 1), "%)")) +
  ylab(paste0("PC2 (", round(pc[2] * 100, 1), "%)")) +
  theme(legend.position = "none", 
        panel.background = element_blank(),
        panel.grid = element_blank()) +
  scale_color_manual(values = c("black", "dark gray"))

p_ploidy

ggsave(paste0("~/cg_andropogon/Figures/traits_pca.byPloidy.", Sys.Date(), ".jpeg"),
       units = "in",
       width = 4,
       height = 4)

# ggarrange(p_group, p_ploidy, ncol = 2, labels = c("A", "B"), font.label = list(size = 22) )
# ggsave(paste0("~/cg_andropogon/Figures/traits_pca.byPloidyandGroup.", Sys.Date(), ".jpeg"),
#        units = "in",
#        width = 14,
#        height = 7)
