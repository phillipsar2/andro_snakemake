### Title: Big bluestem kinship matrix
### Author: Alyssa Phillips
### Date: 1/3/2023

# Code adapted from script provided by Julianna Porter and Jeffrey Ross-Ibarra
# Create a kinship matrix using eqn 4.9 from Caballero (pg 74) - Van Raden (2008) approach

library(ggplot2)
library("pheatmap")
library(RColorBrewer)
library(stringr)
library(viridis)
library(ggbreak) 
library(dplyr)

# 1. Load data ----

##  Genotypes ----

## Common garden files
# geno_tab <- read.table(gzfile("~/Andropogon/kinship/cg.andro.lowcov.nomiss.ibs.gz"), header= T) # no missing data
# geno_tab <- read.table(gzfile("~/Andropogon/pca/lowcov/cg.lowcov.50k.ibs.gz"), header= T) # 50k SNPs

## All Andropogon files
geno_tab <- read.table("~/Andropogon/kinship/all.andro.lowcov.100k.miss20.ibs", header = T)

head(geno_tab)
dim(geno_tab)
### rows = loci
### col = ind

## Load metadata ----
# meta <- read.csv("~/Andropogon/pca/lowcov/cg.geno_meta_02152023.csv", header = T) # common garden
meta <- read.csv("~/Andropogon/pca/lowcov/geno_meta_09212023.csv", header = T) # all andro
head(meta)
dim(meta)
str(meta)

meta$ploidy <- as.factor(meta$ploidy)
meta$population <- as.factor(meta$population)
meta$pop_ploidy <- paste0(meta$population, "_", meta$ploidy)

# 2. Calculate kinship matrix ----

## Convert dataframe to matrix ----
og_mat <- as.matrix(geno_tab[,5:dim(geno_tab)[2]])

## Set missing data to NA ----
paste0("sites with missing data: ", sum(og_mat == -1))

og_mat[og_mat == -1] <- NA

## (1) Calculate snp means accounting for missing data (qs) ----
'qs = average frequency of marker s for all individiuals (N) across populations'
snp_means <- rowSums(og_mat, na.rm = T) / ( dim(og_mat)[2] - rowSums(is.na(og_mat)))
hist(snp_means)

## (2) Replace NAs with the SNP means ----
n = dim(og_mat)[2] # number of genotypes

miss_geno <- colSums(is.na(og_mat)) / dim(og_mat)[1]
hist(miss_geno)

# replace
for (i in 1:dim(og_mat)[1]){
  og_mat[i, is.na(og_mat[i,])] <- snp_means[i]
}

sum(is.na(og_mat))
og_mat[1:10,1:5]

## (3) Calculate sum of squares (denominator of fMij) ----
ssq = sum(snp_means * (1 - snp_means)) 

## (4) Calculate molecular coancestry (fMij) between individuals i and j  ---- 
kin <- matrix(nrow =  n, ncol = n)

for (ind_i in 1:n){
  for (ind_j in 1:n){
    f <- sum( (og_mat[,ind_i] - snp_means) * (og_mat[,ind_j] - snp_means) ) / ssq # ZZ' / ssq
    kin[ind_i, ind_j] <- f
  }
}

kin[1:5,1:5]

## (5) Set diagonal ----

### (A) Center mean on 1 ----
# diag_name = "centereddiag"
# fin_kin <-  kin / mean(diag(kin))
# dim(fin_kin)
# fin_kin[1:5,1:5]

### (B) Set diagonal to NA ----
diag_name = "NAdiag"
fin_kin <- kin
diag(fin_kin) <- NA

### (C) Load in Fis estimates from ngsF ----
## Convert to matrix
# run1_mat <- as.matrix(run1[,5:dim(run1)[2]])
# run2_mat <- as.matrix(run2[,5:dim(run2)[2]])
# colnames(run1_mat) <- r1_name
# colnames(run2_mat) <- r2_name
# 
# # Replace missing data value with NA
# run1_mat[run1_mat == -1] <- NA
# run2_mat[run2_mat == -1] <- NA
# 
# dim(run1_mat)
# dim(run2_mat)
# 
# # Sort so columns are in the same order
# run1_sort <- run1_mat[ ,order(colnames(run1_mat))] 
# run2_sort <- run2_mat[ ,order(colnames(run2_mat))] 
# 
# # Only keep sites that have data in both matrices
# run1_sort[is.na(run2_sort)] <- NA
# run2_sort[is.na(run1_sort)] <- NA
# 
# sum(is.na(run1_sort))
# sum(is.na(run2_sort)) # should have same # of NAs
# 
# miss <- colSums(is.na(run1_sort)) / dim(run1_sort)[1]
# hist(miss*100, 
#      main = "Missing data per genotype") # some genotypes have 80% missing data

## Count the number of sites that drew the same allele in each individual
# # Same read = 1; different reads = 0, no data = NA
# same_read <- ifelse(run1_sort == run2_sort, yes = 1, no = 0)
# 
# hom <- colSums(same_read, na.rm = T) / colSums(!is.na(same_read))
# length(hom)
# 
# hist(hom)
# plot(miss, hom)
# cor(miss, hom)

# fin_kin <- kin
# diag(fin_kin) <- inbred

## missing data versus inbreeding coef
# plot(colSums(is.na(in_af)) / 100000, inbred,
     # xlab = "Missing data per genotype",
     # ylab = "inbreeding coef")
# cor(colSums(is.na(in_af)) / 100000, inbred)

# Is mean associated with missing data?
# plot(miss_geno, diag(fin_kin))
# cor(miss_geno, diag(fin_kin))

## (6) Set values < 0 to 0 ----
fin_kin[fin_kin < 0] <- 0

# 3. Clean up kinship matrix ----

# Add genotypes to columns and rows
colnames(fin_kin) = as.factor(meta$pop_ploidy)
rownames(fin_kin) = as.factor(meta$pop_ploidy)
# colnames(fin_kin) = as.factor(meta$geno_short)
# rownames(fin_kin) = as.factor(meta$geno_short)

# 4. Identify clones ----
colnames(fin_kin) = as.factor(meta$geno_short)
rownames(fin_kin) = as.factor(meta$geno_short)

cutoff = 0.4
sum(fin_kin >= cutoff, na.rm = T) # number of comparisons between clones

potential_clones = fin_kin
potential_clones[!(fin_kin >= cutoff)] <- NA
# View(potential_clones)

only_clones <- apply(potential_clones, 1, function(x) all(is.na(x)))
clones <- potential_clones[ !only_clones, ]
rownames(clones) # These are our clones 

clones_to_drop = c("Afton_TX_3", "Afton_TX_1")
clones_to_keep = c("Afton_TX_2")

meta$genotype[meta$geno_short %in% clones_to_drop]
# [1] "AN20NSCR000409_CKDL210007203-1a-AK4985-AK5009_H3MTWDSX2_L1"    
# [2] "INCU_PCRfree_Texas-1_40_ATCCACTG_Andropogon_gerardii.subsample"

## Drop clones ----
meta_noclones <- meta[!meta$geno_short %in% clones_to_drop,]
dim(meta_noclones)



# 5. Export kinship matrix ----
# write.table(fin_kin, paste0("~/Andropogon/kinship/cg.50k.kinship_mat.", diag_name,".",Sys.Date(),".txt"), row.names = T, col.names = T, quote = F)

# write.table(fin_kin_noclones, paste0("~/Andropogon/kinship/all.andro.lowcov.100k.miss20.", diag_name,".",Sys.Date(),".noclones.txt"), row.names = T, col.names = T, quote = F)

# 6. Plot ----
## Read in exported file ----
# k <- read.table("~/Andropogon/kinship/all.andro.lowcov.100k.miss20.centereddiag.2023-09-27.txt") # all andor
fin_kin_noclones <- read.table("~/Andropogon/kinship/all.andro.lowcov.100k.miss20.NAdiag.2024-03-07.noclones.txt") # no clones for plotting

fin_kin_noclones <- fin_kin[colnames(fin_kin) %in% meta_noclones$geno_short, rownames(fin_kin) %in%  meta_noclones$geno_short] ## Drop clones
dim(fin_kin_noclones)

colnames(fin_kin_noclones) <- rownames(fin_kin_noclones)
fin_kin_noclones <- as.matrix(fin_kin_noclones)
fin_kin_noclones[1:5,1:5]
dim(fin_kin_noclones)

## (1) All samples with clones----
# pdf(file = paste0("~/Andropogon/kinship/cg.50k.kinship_mat.",Sys.Date(),".pdf"),
#     width = 10,
#     height = 9)

# pos_df = data.frame("Pos", meta$population)
p <- pheatmap(fin_kin,  #fin_kin
         fontsize = 6
         )
p
# dev.off()

# ggsave(paste0("~/Andropogon/kinship/cg.50k.kinship_mat.",Sys.Date(),".jpeg"), plot = p )
# ggsave(paste0("~/Andropogon/kinship/all.andro.100k.kinship_mat.",Sys.Date(),".jpeg"),
#        plot = p,
#        height = 15,
#        width = 15,
#        unit = "in")

## (2) All samples without clones ----
colnames(fin_kin_noclones) = as.factor(meta_noclones$geno_short)
rownames(fin_kin_noclones) = as.factor(meta_noclones$geno_short)

mini_pop <- data.frame(
  population = meta_noclones$pop_short
  # population = meta_noclones$genotype
  # cytotype = as.factor(meta_noclones$realploidy)
              )
rownames(mini_pop) <- as.factor(meta_noclones$geno_short)

mini_cyt <- data.frame(
  # population = meta_noclones$pop_short,
  Cytotype = as.factor(meta_noclones$realploidy)
)
rownames(mini_cyt) <- as.factor(meta_noclones$geno_short)

my_colour = list(
  Cytotype = c('6' = "lightgray", '9' = "#707070")
  # random = c(random1 = "#82ed82", random2 = "#9e82ed"),
  # cluster = c(cluster1 = "#e89829", cluster2 = "#cc4ee0")
)

breaks = seq(0,0.5,0.02)

p <- pheatmap::pheatmap(fin_kin_noclones, # fin_kin_noclones
              # fontsize = 12,
              fontsize_row = 8,
              fontsize_col = 8,
              annotation_col = mini_cyt,
              annotation_colors = my_colour,
              annotation_legend = F,
              color = rev(mako(length(breaks))),
              # show_rownames = F,
              # show_colnames = F,
              labels_col = mini_pop$population,
              labels_row = mini_pop$population,
              na_col = "white")

p


# ggsave(paste0("~/Andropogon/kinship/all.andro.100k.kinship_mat.",Sys.Date(),".noclones.jpeg"),
#        plot = p,
#        height = 15,
#        width = 15,
#        unit = "in")

## Subset to KON, BOU, AUS only ----
mix_pops <- c("KON", "BOU", "AUS")
mix_genos <- meta_noclones$geno_short[meta_noclones$pop_short %in% mix_pops]
fin_kin_mix <- fin_kin_noclones[colnames(fin_kin_noclones) %in% mix_genos, rownames(fin_kin_noclones) %in% mix_genos]
dim(fin_kin_mix)

mini_pop <- data.frame(Population = meta_noclones$pop_short[meta_noclones$pop_short %in% mix_pops])
rownames(mini_pop) <- as.factor(meta_noclones$geno_short[meta_noclones$pop_short %in% mix_pops])

mini_cyt <- data.frame(Cytotype = as.factor(meta_noclones$realploidy[meta_noclones$pop_short %in% mix_pops]))
rownames(mini_cyt) <- as.factor(meta_noclones$geno_short[meta_noclones$pop_short %in% mix_pops])
my_colour = list(
  Cytotype = c('6' = "lightgray", '9' = "#707070")
)

breaks = seq(0,0.5,0.02)

mix_p <- pheatmap::pheatmap(fin_kin_mix, # fin_kin_noclones
                        # fontsize = 12,
                        fontsize_row = 10,
                        fontsize_col = 10,
                        color = rev(mako(length(breaks))),
                        annotation_col = mini_cyt,
                        annotation_colors = my_colour,
                        # annotation_row = mini_pop,
                        annotation_legend = T,
                        breaks = breaks,
                        legend_labels  = breaks,
                        # show_rownames = F,
                        # show_colnames = F,
                        labels_col = mini_pop$Population,
                        labels_row = mini_pop$Population,
                        na_col = "white")

ggsave(paste0("~/Andropogon/kinship/mixedpopsonly.andro.100k.kinship_mat.",Sys.Date(),".noclones.jpeg"),
       plot = mix_p,
       height = 7,
       width = 9,
       unit = "in")

## Histograms
library(matrixcalc)
library(reshape2)
kin_noclones_df <- melt(matrixcalc::upper.triangle(fin_kin_noclones))
dim(kin_noclones_df)
head(kin_noclones_df)

nine <- meta_noclones$geno_short[meta_noclones$ploidy == 3]

kin_noclones_df$ploidy1 <- ifelse(kin_noclones_df$Var1 %in% nine, "9x", "6x")
kin_noclones_df$ploidy2 <- ifelse(kin_noclones_df$Var2 %in% nine, "9x", "6x")
head(kin_noclones_df)
kin_noclones_df$combo <- paste0(kin_noclones_df$ploidy1, "_", kin_noclones_df$ploidy2)
kin_noclones_df$combo[kin_noclones_df$combo == "9x_6x"] <- "6x_9x"

hist(kin_noclones_df$value[kin_noclones_df$combo == "9x_9x"])
hist(kin_noclones_df$value[kin_noclones_df$combo == "6x_6x"])
hist(kin_noclones_df$value[kin_noclones_df$combo == "6x_9x"])

hist( sqrt(kin_noclones_df$value) )
hist(kin_noclones_df$value )

ggplot(kin_noclones_df, aes(x = factor(combo), y = sqrt(value))) +
  geom_violin(draw_quantiles = c(.25, .5, .75, .95)) 


kin_noclones_df %>%
  dplyr::filter(combo == "6x_9x" ) %>%
ggplot(aes(fill = factor(combo), x = value)) +
  geom_histogram(fill = "white", color = "black") +
  theme_bw() +
  ylab("Count") + xlab("Kinship") +
  ggtitle("6x-9x")
  # facet_wrap(~combo)
  # scale_x_break(c(0.1, 0.15), scales = 1)
  # facet_wrap(vars(combo))

str(kin_noclones_df)

# install_github("jokergoo/ComplexHeatmap")
# library(ComplexHeatmap)
# cp <- ComplexHeatmap::pheatmap(k, # fin_kin_noclones
#               fontsize = 12,
#               fontsize_row = 7,
#               fontsize_col = 7,
#               annotation_col = mini_cyt,
#               annotation_colors = my_colour,
#               # show_rownames = F,
#               # show_colnames = F,
#               labels_col = mini_pop$population,
#               labels_row = mini_pop$population,
#               na_col = "gray",
#               heatmap_legend_param = list(title = c("Kinship"),
#                                           legend_height = unit(3, "in")
#                                           ),
#               
#               )

# cp





## (3) West samples without 6x clones ----
# Set names as pop_ploidy
colnames(fin_kin_noclones) = as.factor(meta_noclones$geno_short)
rownames(fin_kin_noclones) = as.factor(meta_noclones$geno_short)

# Subset samples to West
# west <- c("Boulder_CO_3","Boulder_CO_2","ESL_2","Konza_3","Konza_2","NM_2","Austin_TX_3", "SAL_2","CDB_2","WEB_2","REL_2","Austin_TX_2","BB_2","Afton_TX_2", "NE_2","SAL_3")
# meta_west <- meta_noclones[meta_noclones$pop_ploidy %in% west,]
# west_genos <- meta_west$geno_short
# length(west_genos)

pop_assignment <- read.csv("~/Andropogon/structure/no_clones/pop_assignment2024-03-27.csv", header = T)
west <- pop_assignment$geno[pop_assignment$group == "1"]
meta_west <- meta_noclones[meta_noclones$geno_short %in% west,]
west_genos <- meta_west$geno_short
length(west_genos)


# fin_kin_sub <- fin_kin_noclones[colnames(fin_kin_noclones) %in% west,rownames(fin_kin_noclones) %in% west]
fin_kin_sub <- fin_kin_noclones[colnames(fin_kin_noclones) %in% west_genos, rownames(fin_kin_noclones) %in% west_genos]
dim(fin_kin_sub)

# Plot
mini_pop_W <- data.frame(population = meta_west$pop_short)
rownames(mini_pop_W) <- as.factor(meta_west$geno_short)

mini_cyt_W <- data.frame(Cytotype = as.factor(meta_west$realploidy))
rownames(mini_cyt_W) <- as.factor(meta_west$geno_short)

my_colour_W = list(Cytotype = c('6' = "lightgray", '9' = "#707070"))

breaks = seq(0,0.5,0.02)

p <- pheatmap::pheatmap(fin_kin_sub, 
              fontsize = 12,
              annotation_col = mini_cyt_W,
              annotation_colors = my_colour_W,
              annotation_legend = F,
              color = rev(mako(length(breaks))),
              labels_col = mini_pop_W$population,
              labels_row = mini_pop_W$population,
              legend_labels  = breaks,
              breaks = breaks,
              na_col = "white")
p


ggsave(paste0("~/Andropogon/kinship/west.andro.100k.kinship_mat.", Sys.Date(),".noclones.jpeg"),
       plot = p,
       height = 15,
       width = 15,
       unit = "in")

# pl_color <- as.matrix(x = NA, nrow = length(meta$ploidy), ncol = 1)
# pop_color <- as.matrix(x = NA, nrow = length(meta$population), ncol = 1)
# 
# for (i in 1:length(meta$ploidy)){
#   if (meta$ploidy[i] == 2){
#     pl_color[i] <- "white"
#   }
#   else{
#     pl_color[i] <- "black"
#   }
# }
# 
# 
# heatmap(fin_kin,
#         RowSideColors = pl_color,
#         # ColSideColors = pop_color,
#         xlab = "population"
# )

# 7. Calculate statistics ----
# Overall mean
fin_kin_noclones[upper.tri(fin_kin_noclones, diag = F)] %>% mean()
fin_kin_noclones[upper.tri(fin_kin_noclones, diag = F)] %>% sd()

colnames(fin_kin_noclones) <- rownames(fin_kin_noclones)
ind <- which( upper.tri(fin_kin_noclones,diag=F) , arr.ind = TRUE )

kin_noclones_df <- data.frame( geno1 = dimnames(fin_kin_noclones)[[2]][ind[,2]] ,
            geno2 = dimnames(fin_kin_noclones)[[1]][ind[,1]] ,
            kinship = fin_kin_noclones[ ind ] )
dim(kin_noclones_df)
head(kin_noclones_df)

# Among 9x
geno9 <- meta$geno_short[meta_noclones$ploidy == "3"]

kin_9x <- filter(kin_noclones_df, geno1 %in% geno9 & geno2 %in% geno9)
mean(kin_9x$kinship, na.rm = T)
sd(kin_9x$kinship, na.rm = T)

max(kin_9x$kinship)

# Among 6x
geno6 <- meta$geno_short[meta_noclones$ploidy == "2"]

kin_6x <- filter(kin_noclones_df, geno1 %in% geno6 & geno2 %in% geno6)
mean(kin_6x$kinship, na.rm = T)
sd(kin_6x$kinship, na.rm = T)

# Between 6x and 9x

