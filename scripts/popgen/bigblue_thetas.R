### Title: Calculating thetas of the 6x low-coverage samples
### Author: Alyssa Phillips
### Date: 7/18/2023

library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)
library(gridExtra)
library(paletteer)
library(gtable)
library(stringr)
library(ggpubr)

pop_names <- str_split(list.files("~/Andropogon/thetas/raw_thetas/"), "[.]", simplify = T)[,1] %>% unique()
window_size <- str_split(list.files("~/Andropogon/thetas/raw_thetas/"), "[.]", simplify = T)[,3] %>% unique()

# Calculate per-site thetas ----
# pop = population ID (character)
# window = window size (numeric)
# theta = c("tW", "tP", "Tajima")
calc_persite <- function(pop, window, theta){
  df <- read.csv(paste0("~/Andropogon/thetas/raw_thetas/", pop, ".all.", window, ".thetas.idx.pestPG.gz"), sep = "\t")
  sub_df <- df %>% filter(df$nSites / window >= 0.1) # remove windows with < 10% of sites
  if(theta == "tP"){ # calculate per-site theta in each window 
    t <- sub_df$tP / sub_df$nSites 
  } else if(theta == "tW"){
    t <- sub_df$tW / sub_df$nSites
  } else if(theta == "Tajima"){
    # t <- sub_df$Tajima / sub_df$nSites # divide by sites
    t <- sub_df$Tajima # or don't divide by sites
  }
  global_t <- sum(t)/dim(sub_df)[1] # mean theta
  sd_t <- sd(t) # standard deviation of theta
  c(pop, global_t, (global_t - sd_t), (global_t + sd_t)) # output
}

## tP ----
tp_10 <- matrix(nrow = length(pop_names), ncol = 4)
colnames(tp_10) <- c("pop", "tP", "minsd", "maxsd")
for (i in 1:length(pop_names)){
  tp_10[i,] <- calc_persite(pop_names[i], 10000, "tP")
} # all different length

tp_50 <- matrix(nrow = length(pop_names), ncol = 4)
colnames(tp_50) <- c("pop", "tP", "minsd", "maxsd")
for (i in 1:length(pop_names)){
  tp_50[i,] <- calc_persite(pop_names[i], 50000, "tP")
} # all different length

## tW ----
tw_10 <- matrix(nrow = length(pop_names), ncol = 4)
colnames(tw_10) <- c("pop", "tW", "minsd", "maxsd")
for (i in 1:length(pop_names)){
  tw_10[i,] <- calc_persite(pop_names[i], 10000, "tW")
} # all different length

tw_50 <- matrix(nrow = length(pop_names), ncol = 4)
colnames(tw_50) <- c("pop", "tW", "minsd", "maxsd")
for (i in 1:length(pop_names)){
  tw_50[i,] <- calc_persite(pop_names[i], 50000, "tW")
} # all different length

## Tajima's D ----
# Consideration on autopolyploidy: https://www.frontiersin.org/articles/10.3389/fgene.2018.00480/full

td_10 <- matrix(nrow = length(pop_names), ncol = 4)
colnames(td_10) <- c("pop", "Tajima", "minsd", "maxsd")
for (i in 1:length(pop_names)){
  td_10[i,] <- calc_persite(pop_names[i], 10000, "Tajima")
} # all different length

td_50 <- matrix(nrow = length(pop_names), ncol = 4)
colnames(td_50) <- c("pop", "Tajima", "minsd", "maxsd")
for (i in 1:length(pop_names)){
  td_50[i,] <- calc_persite(pop_names[i], 50000, "Tajima")
} # all different length


# Plot ----
# scaleFUN <- function(x) sprintf("%.2f", x)

## Average tP, 10k windows ----
tp_10_df <- as.data.frame(tp_10)
tp_10_df[,2:4] <- lapply(tp_10_df[,2:4], as.numeric)

# A
p_tp10 <- tp_10_df %>%
  ggplot() +
  geom_segment( aes(x = pop, xend = pop, y=minsd, yend=maxsd), color="grey") +
  geom_point( aes(x = pop, y = tP), color = "black", size=2 ) +
  # geom_point( aes(x = pop, y = minsd), color = "gray", size=2 ) +
  # geom_point( aes(x = pop, y = maxsd), color = "gray", size=2 ) +
  coord_flip()+
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank()
    # panel.border = element_blank(),
  ) +
  xlab("Population") +
  ylab("Nucleotide diversity") +
  ylim(0.003, 0.01)
  # scale_y_continuous(labels=scaleFUN)

## Average tP, 50k windows ----
tp_50_df <- as.data.frame(tp_50)
tp_50_df[,2:4] <- lapply(tp_50_df[,2:4], as.numeric)

# D
p_tp50 <- tp_50_df %>%
  ggplot() +
  geom_segment( aes(x = pop, xend = pop, y=minsd, yend=maxsd), color="grey") +
  geom_point( aes(x = pop, y = tP), color = "black", size=2 ) +
  # geom_point( aes(x = pop, y = minsd), color = "gray", size=2 ) +
  # geom_point( aes(x = pop, y = maxsd), color = "gray", size=2 ) +
  coord_flip()+
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    # axis.title.y = element_blank()
    # panel.border = element_blank(),
  ) +
  xlab("Population") +
  ylab("Nucleotide diversity") + 
  ylim(0.003, 0.01)
# scale_y_continuous(labels=scaleFUN)

p_tp <- ggarrange(p_tp10, p_tp50, ncol = 1, labels = c('A', 'D'))

# ggsave(paste0("~/Andropogon/thetas/plots/nucleotide_div.6x_lowcov.",Sys.Date(),".jpeg"), 
#        plot = p_tp, 
#        dpi = 300, 
#        height = 6, 
#        width = 6, unit = "in")

## Average tW, 10k windows ----
tw_10_df <- as.data.frame(tw_10)
tw_10_df[,2:4] <- lapply(tw_10_df[,2:4], as.numeric)

# B
p_tw10 <- tw_10_df %>%
  ggplot() +
  geom_segment( aes(x = pop, xend = pop, y=minsd, yend=maxsd), color="grey") +
  geom_point( aes(x = pop, y = tW), color = "black", size=2 ) +
  # geom_point( aes(x = pop, y = minsd), color = "gray", size=2 ) +
  # geom_point( aes(x = pop, y = maxsd), color = "gray", size=2 ) +
  coord_flip()+
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.title.x = element_blank()
    # panel.border = element_blank(),
  ) +
  xlab("Population") +
  ylab("Watterson's theta") +
  # ylim(0.003, 0.011) +
  scale_y_continuous(breaks = seq(0.004,0.01, 0.002))
# scale_y_continuous(labels=scaleFUN)

## Average tW, 50k windows ----
tw_50_df <- as.data.frame(tw_50)
tw_50_df[,2:4] <- lapply(tw_50_df[,2:4], as.numeric)

# E
p_tw50 <- tw_50_df %>%
  ggplot() +
  geom_segment( aes(x = pop, xend = pop, y=minsd, yend=maxsd), color="grey") +
  geom_point( aes(x = pop, y = tW), color = "black", size=2 ) +
  # geom_point( aes(x = pop, y = minsd), color = "gray", size=2 ) +
  # geom_point( aes(x = pop, y = maxsd), color = "gray", size=2 ) +
  coord_flip()+
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank()
    # panel.border = element_blank(),
  ) +
  xlab("Population") +
  ylab("Watterson's theta") +
  # ylim(0.003, 0.011) +
  scale_y_continuous(breaks = seq(0.004,0.01, 0.002))

p_tw <- ggarrange(p_tw10, p_tw50, ncol = 1, labels = c('B', 'E'), hjust = 0.5)

# ggsave(paste0("~/Andropogon/thetas/plots/wattersons_theta.6x_lowcov.",Sys.Date(),".jpeg"), 
#        plot = p_tw, 
#        dpi = 300, 
#        height = 6, 
#        width = 6, unit = "in")

## Average Tajima's D, 10k windows ----
td_10_df <- as.data.frame(td_10)
td_10_df[,2:4] <- lapply(td_10_df[,2:4], as.numeric)

# C
p_td10 <- td_10_df %>%
  ggplot() +
  geom_segment( aes(x = pop, xend = pop, y=minsd, yend=maxsd), color="grey") +
  geom_point( aes(x = pop, y = Tajima), color = "black", size=2 ) +
  # geom_point( aes(x = pop, y = minsd), color = "gray", size=2 ) +
  # geom_point( aes(x = pop, y = maxsd), color = "gray", size=2 ) +
  coord_flip()+
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.title.x = element_blank()
    # panel.border = element_blank(),
  ) +
  xlab("Population") +
  ylab("Tajima's D (10 Kbp windows)") +
  scale_y_continuous(limits = c(-2, 2))       
  # scale_y_continuous(labels = scales::label_number(), limits = c(-0.0008, 0.0012))


## Average Tajima's D, 50k windows ----
td_50_df <- as.data.frame(td_50)
td_50_df[,2:4] <- lapply(td_50_df[,2:4], as.numeric)

# F
p_td50 <- td_50_df %>%
  ggplot() +
  geom_segment(aes(x = pop, xend = pop, y=minsd, yend=maxsd), color="grey", size = 0.5) +
  geom_point( aes(x = pop, y = Tajima), color = "black", size=2 ) +
  # geom_point( aes(x = pop, y = minsd), color = "gray", size=2 ) +
  # geom_point( aes(x = pop, y = maxsd), color = "gray", size=2 ) +
  coord_flip()+
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none", 
    axis.title.y = element_blank()
    # panel.border = element_blank(),
  ) +
  xlab("Population") +
  ylab("Tajima's D") +
  ylim(-2, 2) 
  # scale_y_continuous(labels = scales::label_number(),
  #                    limits = c(-0.0008, 0.0012)) 

p_td <- ggarrange(p_td10, p_td50, 
                  ncol = 1, 
                  labels = c('C', 'F'),
                  hjust = 0.5
                  )

# ggsave(paste0("~/Andropogon/thetas/plots/tajimas_D.6x_lowcov.per_window.",Sys.Date(),".jpeg"),
#        plot = p_td,
#        dpi = 300,
#        height = 6,
#        width = 6, unit = "in")


## Combined image ----

p_all <- ggarrange(p_tp, p_tw, p_td, ncol = 3)

ggsave(paste0("~/Andropogon/thetas/plots/thetas.6x_lowcov.",Sys.Date(),".jpeg"),
       plot = p_all,
       dpi = 300,
       height = 5,
       width = 12, unit = "in")

# Species-wide averages table ----
avgs <- cbind(td_10_df$Tajima, td_50_df$Tajima, tp_10_df$tP, tp_50_df$tP, tw_10_df$tW, tw_50_df$tW) %>% colMeans() %>% as.data.frame()
rownames(avgs) <- c("Tajima10", "Tajima50", "tP10", "tP50", "tW10", "tW50")
avgs

########## Old scripts #################
# # Add pi from other species
# species <- c("maize", "A. thaliana", "loblolly pine")
# lower_pi <- c(.006, .0011, 0.00027)
# upper_pi <- c(.01,.002, 0.01728)
# 
# other <- as.data.frame(cbind(species, lower_pi, upper_pi))
# str(other)
# other$upper_pi <- as.numeric(as.character(other$upper_pi))
# other$lower_pi <- as.numeric(as.character(other$lower_pi))
# 
# g2 <- ggplot(other, aes(color=species)) +
#   geom_segment( aes(x=species, xend=species, y=lower_pi, yend=upper_pi)) +
#   geom_point( aes(x=species, y=lower_pi), size=3) +
#   geom_point( aes(x=species, y=upper_pi), size=3 ) +
#   coord_flip()+
#   theme_ipsum(axis_title_just = "cc" ) + 
#   scale_color_paletteer_d("nationalparkcolors::CraterLake", -1) +
#   theme(legend.position = "none", 
#         axis.title.y = element_blank(),
#         plot.margin = unit(c(0,10,10,10), "pt")) +
#   ylab("per-site nucleotide diversity") +
#   ylim(0,0.04)
# 
# 
# p1 <- ggplotGrob(g1)
# p2 <- ggplotGrob(g2)
# p <- rbind(p1,p2, size = "first")
# #library(grid)
# p$widths <- unit.pmax(p1$widths, p2$widths)
# 
# grid.newpage()
# grid.draw(p)
# 
# ggsave("~/poa/pi/multiple_species_pi_09142021.jpg",plot = p, dpi = 300, height = 15, width = 25, unit = "cm")
