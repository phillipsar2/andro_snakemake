### Title: STRUCTURE
### Author: Alyssa Phillips
### Date: 1/12/2023

###
### Create STRUCTURE input ----
### 

library(dplyr)

# (1) Load data ----
## Common garden
# g <- read.table(gzfile("~/Andropogon/pca/lowcov/cg.lowcov.50k.ibs.gz"), header= T) # 50k SNPs

## All andro
g <- read.table("~/Andropogon/kinship/all.andro.lowcov.100k.miss20.ibs", header = T) # 100k SNPs
head(g)
dim(g)
str(g)

# Load meta data
## Common garden
# meta <- read.csv("~/Andropogon/pca/lowcov/cg.geno_meta_01102023.csv", header = T)
## All Andro
meta <- read.csv("~/Andropogon/pca/lowcov/geno_meta_09212023.csv")
head(meta)
dim(meta)
str(meta)

# (2) Prep data ----
## Exclude clones ----
clones_to_drop = c("Afton_TX_3", "Afton_TX_1")
clones_to_drop_index = meta$index[meta$geno_short %in% clones_to_drop]

meta_noclones <- meta[!meta$geno_short %in% clones_to_drop,]
dim(meta_noclones)

g_noclones <- g[,!colnames(g) %in% clones_to_drop_index]
dim(g_noclones)

## Replace -1 with -9 (missing data marker) ----
g_noclones[g_noclones == -1] <- -9
sum(g_noclones == -9) # missing data

## Assign each pop and ID ----
npops <- length(unique(meta_noclones$pop_short))
index <- seq(from = 1, to = npops, by = 1)
pop_index <- data.frame(pop = unique(meta_noclones$pop_short), index = index)

full_pop_index <- pop_index$index[match(x = meta_noclones$pop_short, table = pop_index$pop)]

# (3) Create table ----
out <- data.frame(meta_noclones$geno_short, # optional genotype info
                  full_pop_index, # optional population index info
                  t(g_noclones[,5:dim(g_noclones)[2]])) # required genotype data
dim(out)
out[1:5,1:5]

# (4) Export table ----
write.table(out, file = paste0("~/Andropogon/structure/all.andro.noclones.lowcov.100k.", Sys.Date(),".structure_input.txt"),
            sep = " ",
            row.names = F,
            col.names = F,
            quote = F)

write.table(pop_index, file = paste0("~/Andropogon/structure/all.andro.noclones.lowcov.100k.",Sys.Date(),".pop_index.txt"), 
            quote = F, 
            row.names = F)


###
### Plot STRUCTURE output ----
###

library(ggplot2)
library(forcats)
library(ggthemes)
library(patchwork)
library(tidyr)
library(dplyr)

# Turn structure output into a new file 
k = 2
run = 1

## Common garden 
# pop_index <- read.table("~/Andropogon/structure/cg.lowcov.50k.pop_index.txt", header = T)
# str <- read.table(paste0("~/Andropogon/structure/common_garden/ancestry_only/cg.lowcov.50k.k",k,".run",run,".75steps.ancestry.txt"), 
#                    header = T)
# str[str$geno == "INCR_PCRfre", 2:3 ] <- c("WAL_3", "12") ## replace the mislabeled 12mi_3 with WAL_3

## All Andro - without clones
pop_index <- read.table("~/Andropogon/structure/all.andro.noclones.lowcov.100k.2024-03-08.pop_index.txt", header = T)
str <- read.table(paste0("~/Andropogon/structure/no_clones/ancestry/all.andro.noclones.lowcov.100k.k",k,".run",run,".75steps.ancestry.txt"),
                  header = T)

## All Andro except inbreds
# pop_index <- read.table("~/Andropogon/structure/all.andro.lowcov.100k.pop_index.txt", header = T)
# str <- read.table(paste0("~/Andropogon/structure/all_andro/ancestry_only/noinbreds.andro.lowcov.100k.k",k,".run",run,".ancestry.txt"), 
#                   header = T, sep = " ", row.names = NULL)

colnames(str) <- c("index", "geno", "Miss", "assigned_pop", paste0("q", seq(1, k)))
head(str)

## Convert wide to long
str_long <- gather(str, cluster, q, paste0("q", seq(1, k)), factor_key = T)
str_long$index <- as.factor(str_long$index)
str_long$assigned_pop <- as.factor(str_long$assigned_pop)

str_long$pop <- NA
for (i in 1:dim(str_long)[1]){
  str_long$pop[i] <- pop_index[pop_index$index == str_long$assigned_pop[i],]$pop
}
head(str_long)
tail(str_long)

# str_long <- str_long[,-5]
# head(str_long)
# str_long$pop <- pop_index[match(str_long$assigned_pop, pop_index$index),1]

# > Export group assignment ----
pop_assignment <- c()
for (i in 1:dim(str)[1]){
  pop_assignment[i] <- ifelse(str$q1[i] > str$q2[i], "1", "2")
}

pop_assignment <- cbind(str$geno,pop_assignment )
colnames(pop_assignment) <- c("geno", "group")

# write.csv(pop_assignment,paste0("~/Andropogon/structure/no_clones/pop_assignment", Sys.Date(),".csv"), quote = F, row.names = F, col.names = T)

# > Plot ----

## order genotypes by cluster proportion
tmp_lvls <- arrange(str_long, cluster, -q)$index %>% unique()
str_long$index <- factor(str_long$index, levels = tmp_lvls)

str_long$pop <- as.factor(str_long$pop)
str(str_long)

two_group_colors = c("#E1BE6A", "#40B0A6") # A

# plot
means <- aggregate(q ~ pop + cluster, data = str_long, FUN = mean) %>% 
  filter(cluster == "q2")
str_long$pop <- factor(str_long$pop, levels = means$pop[order(means$q)])

arrange(str_long, cluster, -q) %>% # sort genotypes by admix proportion
  na.omit() %>%
  ggplot(aes(x = index,
             y = q, 
             fill = cluster)) +
    geom_col(color = "gray", linewidth = 0.1) +
    # facet_grid(~fct_inorder(pop), switch = "x", scales = "free", space = "free") +
    facet_grid(~ pop, switch = "x", 
               scales = "free", space = "free") +
    theme_minimal() + 
  labs(x = "Individuals", 
       # title = paste0("K = ", k), 
       y = "Ancestry proportion") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = expand_scale(add = 1)) +
    # scale_fill_gdocs(guide = FALSE) + # for k < 10
  # scale_fill_manual(values = two_group_colors) + # for k = 2
    theme(
      panel.spacing.x = unit(0.1, "lines"),
      axis.text.x = element_blank(),
      panel.grid = element_blank(),
      strip.text.x = element_text(size = 16, angle = 90),
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 16),
      axis.text.y = element_text(margin = margin(r = -1), size = 10)
    ) 

# > Save plot  ----
## Common garden 
# ggsave(paste0("~/Andropogon/structure/common_garden/plots/cg.lowcov.50k.k",k,".run",run,".75steps.ancestry.pdf"),
#        width = 10,
#        height = 4)

## All andropogon - no clones
ggsave(paste0("~/Andropogon/structure/no_clones/plots/all.andro.noclones.lowcov.100k.k",k,".run",run,".main.75steps.ancestry.pdf"),
       width = 24, #26
       height = 3) #4

## All andro, no inbreds
# ggsave(paste0("~/Andropogon/structure/all_andro/plots/noinbreds.andro.lowcov.100k.k",k,".run",run,".ancestry.pdf"),
#        width = 25,
#        height = 4)

# Evaluate K ----
## Common garden
# pr <- read.table("~/Andropogon/structure/common_garden/ancestry_only/prob_of_data.75steps.txt", header = T, row.names = NULL)

## All Andropogon
# pr <- read.table("~/Andropogon/structure/all_andro/ancestry_only/prob_of_data.75steps.txt", header = T, row.names = NULL)
# pr <- read.table("~/Andropogon/structure/all_andro/ancestry_only/noinbreds.prob_of_data.txt", header = T, row.names = NULL)
pr <- read.table("~/Andropogon/structure/no_clones/ancestry/all.andro.noclones.prob_of_data.txt", header = T, row.names = NULL)

# Plot log likelihoods - higher = better
ggplot(pr, aes(x = k, y = Pr, color = as.factor(run), group = run)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  ylab(" Ln P(D)") +
  xlab("K") +
  scale_color_manual(name = "Run", values = c("blue", "black", "gray"))

ggsave("~/Andropogon/structure/no_clones/best_K_noclones.20240404.jpeg",
       width = 6, height = 4, unit = "in")


# Admixture pie charts ----
library(adegenet)
library(poppr)
# library(LEA)
library(reshape2)
library(dplyr)
library(ggplot2)
library(rworldmap)
library(rworldxtra)
library(ggsn)
library(sf)
library(raster)
library(rgeos)
library(maps)
library(maptools)
library(grid)
library(miscTools)
library(stringr)
library(ggpubr)

## Load metadata
pop_meta <- read.csv("~/Andropogon/Maps/04242023/AllAndro_meta.csv", header = T)
head(pop_meta)
colnames(pop_meta)[5] <- "population"

## add population data to q data
str$population <- NA
for (i in 1:length(str$assigned_pop)){
  str$population[i] <- pop_index[str$assigned_pop[i],]$pop
}
dim(str)
dim(pop_meta)

## merge metadata with q data
# map_data <- merge(str, pop_meta, "population", all.x = T, sort = F)

## > Calculate mean admixture proportions for each site ----
q1_mean <- aggregate(q1 ~ population, data = str, mean)
q2_mean <- aggregate(q2 ~ population, data = str, mean)

## Merge with lat long data then melt from wide to long
avg_q <- data.frame(merge(q1_mean, pop_meta, "population"),
           q2 = q2_mean$q2)

avg_q_long <- gather(avg_q, cluster, q, c("q1", "q2"), factor_key = T)

## Define a function to plot pie charts using ggplot for each site
## following https://github.com/Tom-Jenkins/admixture_pie_chart_map_tutorial/blob/master/pie_chart_admixture_map_tutorial.R 
pie_charts = function(admix_df, site, cols){
  # admix_df = dataframe in long format of admixture proportions per site 
  # site = string 
  # cols = vector of colours of length(clusters)
  ggplot(data = subset(admix_df, population == site),
         aes(x = "", y = q, fill = cluster))+
    geom_bar(width = 1, stat = "identity", colour = "black", show.legend = FALSE)+
    coord_polar(theta = "y")+
    scale_fill_manual(values = cols)+
    theme_void()
}

## Test function on one site
pie_charts(avg_q_long, site = "TWE", cols = c("#E1BE6A", "#40B0A6"))

## > Apply function to all sites using for loop ----
pies = list()
for (i in unique(avg_q_long$population)){
  pies[[i]] = pie_charts(admix_df = avg_q_long, site = i, cols = c("#E1BE6A", "#40B0A6")) 
}

## > Make basemap ----

# Set map boundary (xmin, xmax, ymin, ymax)
summary(pop_meta)
boundary = extent(-107, -75, 26, 52)
boundary

# Get map outlines from rworldmap package
map.outline = getMap(resolution = "high")

# Crop to boundary and convert to dataframe
map.outline = crop(map.outline, y = boundary) %>% fortify()

# >> Plot basemap ----
basemap = ggplot() +
  geom_polygon(data = map.outline, aes(x = long, y = lat, group = group), fill = "grey",
               colour="black", size=0.5)+
  coord_quickmap(expand  =F) +
  ggsn::north(map.outline, symbol = 10, scale = 0.06, location = "topleft") + # north arrow
  # ggsn::scalebar(data = map.outline, dist = 200, dist_unit = "km", height = 0.01,
  #                transform = TRUE, model = "WGS84", 
  #                location = "bottomleft", anchor = c(x = -12.5, y = 45),
  #                st.bottom = FALSE, st.size = 4, st.dist = 0.015)+
  xlab("Longitude")+
  ylab("Latitude")+
  theme(
    axis.text = element_text(colour="black", size=12),
    axis.title = element_text(colour="black", size=14),
    panel.background = element_rect(fill="lightsteelblue2"),
    panel.border = element_rect(fill = NA, colour = "black", size = 0.5),
    panel.grid.minor = element_line(colour="grey90", size=0.5),
    panel.grid.major = element_line(colour="grey90", size=0.5),
    legend.text = element_text(size=12),
    legend.title = element_blank(),
    legend.key.size = unit(0.7, "cm"),
    legend.position = "top"
  )

basemap

##  > Add pie charts ----

radius = 1 # Define pie chart sizes

# Convert ggplot pie charts to annotation_custom layers
subsites = unique(avg_q_long$population)
coord.list = list()
for (i in subsites){ # create list of coordinates
  coord.list[[i]] = c(subset(pop_meta, population == i)$Longitude, subset(pop_meta, population == i)$Latitude)
}

pies.ac = list()
for (i in 1:length(subsites)){
  pies.ac[[i]] = annotation_custom(grob = ggplotGrob(pies[[i]]),
                                   xmin = coord.list[[i]][[1]] - radius,
                                   xmax = coord.list[[i]][[1]] + radius,
                                   ymin = coord.list[[i]][[2]] - radius,
                                   ymax = coord.list[[i]][[2]] + radius)
}

# Add pies on
pie.map = basemap + pies.ac
pie.map
# ggsave("~/cg_andropogon/Figures/pie_charts_map.20240404.jpeg", width = 8, height = 10, dpi = 300)
