## Potting nQuire results
## Author: Alyssa Phillips
## Date: 10/28/2021

library(dplyr) # v3.6.2
library(mclust) # v5.4.7
library(stringr)

# Read in results from nQuire ----
data <- read.csv("~/Andropogon/nQuire/nquire_results_denoised.txt", header = T, sep = "\t")

levels(data$ploidy) <- c("6x", "9x", "unknown")

data[is.na(data)] <- "unknown"

head(data)
str(data)

# Normalize max log-likelihoods to the free model ----
normL <- data[,4:6] / data[,3]
head(normL)

# Plot ----
'closer to 1 is better'
# plot(normL, upper.panel = NULL)

data.frame(ploidy = data$ploidy, normL) %>%
  ggplot(aes(x = diploid, y = triploid, color = as.factor(ploidy))) +
  geom_point(size = 2) +
  theme_bw() +
  scale_color_manual(values = c("#648FFF","#DC267F","#FFB000" ),
                     labels = c("6x", "9x", "Unknown")) +
  guides(color=guide_legend(title="Ploidy")) +
  xlab("Diploid") +
  ylab("Triploid")

ggsave("~/Andropogon/Figures/nQuire.20240407.jpeg",
       width = 5,
       height = 4,
       unit = "in")

# plot(normL$dip, normL$tri,
#      xlab = "Diploid",
#      ylab = "Triploid")
# plot(normL$tet, normL$tri)
# plot(normL$dip, normL$tet)

# create scatter plot for each pair of variable
# colnames(normL) <- c("diploid", "triploid", "tetraploid")
# clPairs(normL, 
#         data$ploidy, 
#         lower.panel = NULL,)
# clPairsLegend(0.1, 0.3, 
#               class = clp$class, 
#               col = clp$col, 
#               pch = clp$pch)
