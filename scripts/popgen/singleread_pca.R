### Title: Single read PCA script
### Author: Alyssa Phillips
### Date: 9/13/22

library(stringr)
library(ggplot2)

# (0) Load data ----
## 1. Load covariance matrix ----
# 100k snps without clones, all andro
C <- as.matrix(read.table("~/Andropogon/pca/lowcov/all.andro.lowcov.miss20.100k.covMat"))
# 50k snps for all ando
# C <- as.matrix(read.table("~/Andropogon/pca/lowcov/all.andro.lowcov.50k.covMat")) 
# g <- read.table(gzfile("~/Andropogon/pca/lowcov/all.andro.lowcov.ibs.gz"), row.names = NULL, fill = T, header = T)
# head(g)
# g[1:5,1:5]
# dim(g)

# 50k snps, common garden only
# C <- as.matrix(read.table("~/Andropogon/pca/lowcov/cg.lowcov.50k.covMat")) 

# remove duplicates
# dups <- duplicated(pl$row.names)
# pl[pl$row.names != dups, ] %>% dim()
# pl_df <- pl[!duplicated(pl$row.names), ]
# row.names(pl_df) <- pl_df$row.names
# pl_df[1:3,1:3]

# drop row.names column
# pl_df <- pl_df[1:dim(pl_df)[1], 2:dim(pl_df)[2]]
# dim(pl_df)

## 2. Load metadata ----
## has the index for run with all ind and clones
meta <- read.csv("~/Andropogon/pca/lowcov/geno_meta_09212023.csv", header = T)
head(meta)
dim(meta)
str(meta)

# > Drop clones ----
clones_to_drop = c("Afton_TX_3", "Afton_TX_1")
meta_noclones <- meta[!meta$geno_short %in% clones_to_drop,]

# (1)  Single Read PCA ----
# Compute eigenvalues and corresponding eigenvectors of S
# these are the principal components of S
e <- eigen(C)

# Print proportion of total var explained by the components
pc <- e$values / sum(e$values)
pc[1:10]

# Plot with ggplot
df <- as.data.frame(cbind(e$vectors[,1:10], meta_noclones))
names(df) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", colnames(meta_noclones))
df$ploidy <- as.factor(df$ploidy)
head(df)
str(df)

## Common garden pop order (W to E)
# populationOrder <- c("Boulder_CO","NM", "SAL",  "CDB", "REL", "Austin_TX", "ESL", "MBTGPP", "Konza", "FUL", "12mi", "WI", "WAL", "Suther")
# df$population <- factor(df$population, levels = populationOrder)

## All Andropogon pop order (W to E)

g1 <- df[complete.cases(df[,13:14]),] %>%
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_point(aes(fill = pop_short, color = pop_short, shape = ploidy),
             size = 2) +
  xlab(paste0("PC1 (", round(pc[1]*100,1), "%)")) +
  ylab(paste0("PC2 (", round(pc[2]*100, 1), "%)")) +
  theme_bw(base_size = 12) +
  theme(
    # legend.position="none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) + 
  guides(color = guide_legend(title = "Population"),
         shape = 'none',
         fill = 'none')
  # ggtitle("All Andropogon, 100k SNPs")


g1

# ggsave(paste0("~/Andropogon/pca/all.andro.noclones.100k.", Sys.Date(), ".pop.jpeg"),
#        units = "in",
#        width = 6,
#        height = 4)

# (2) Plot by W & E ----
# pop_assignment <- read.csv("~/Andropogon/structure/no_clones/pop_assignment2024-03-21.csv", header = T)
pop_assi_df <- as.data.frame(pop_assignment)

pop_assign <- merge(df, pop_assi_df, by.x = "geno_short", by.y = "geno" )


pop_assign$groupname <- ifelse(pop_assign$group == '1', "West", "East")
str(pop_assign)

# west <- c("AFT", "MIL", "BAR", "BOU", "CDB", "REL", "SAL", "WEB", "KON", "AUS")
# df$group <- c()
# df$group <- ifelse(df$pop_short %in% west, "West", "East")

pop_assign[complete.cases(pop_assign[,17]),] %>%
  ggplot(aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = groupname, shape = ploidy),
             size = 2, alpha = 0.5) +
  xlab(paste0("PC1 (", round(pc[1]*100,1), "%)")) +
  ylab(paste0("PC2 (", round(pc[2]*100, 1), "%)")) +
  theme_bw(base_size = 12) +
  theme(
    # legend.position="none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  ggtitle("All sampled andropogon") +
  # scale_color_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("#40B0A6", "#E1BE6A")) +
  scale_shape(name = "Ploidy", labels = c("6x", "9x")) + 
  guides(color = guide_legend(title = "Group")) +
  ggtitle(NULL)

ggsave(paste0("~/Andropogon/pca/all.andro.noclones.100k.", Sys.Date(), ".groups.jpeg"),
       units = "in",
       width = 5,
       height = 4)
