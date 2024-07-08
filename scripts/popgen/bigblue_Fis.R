### Title: Plotting individual inbreeding coefficients estimated from ngsF
### Author: Alyssa Phillips
### Date: 6/20/23

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(purrr)
library(matrixStats)
library(emmeans)
library(multcomp)

# Load Fis data ----

dir = "~/Andropogon/inbreeding/output"
files <- list.files(path = dir, pattern = "*.allruns.indF") 
files <- paste0(dir,"/",files)
files

# file_list <- c("~/Andropogon/inbreeding/output/highcov_6x_andro.Chr*A.allruns.indF") 
extract_data_from_file <- function(file) { 
  data <- read.csv(file, header = F)
  data <- data[,1] %>% as.data.frame()
  return(data) 
} 

extracted_data_list <- map_dfc(files, 
                              extract_data_from_file
                              )
length(files) == 30 # num of chromosomes
dim(extracted_data_list) # ind x chrom

chr <- str_split(files, pattern = "\\.", simplify = T )[,2]
colnames(extracted_data_list) <- chr

# Load meta data ----

geno_list <- read.csv("~/Andropogon/inbreeding/output/highcov.6x.bamlist", header = F)
colnames(geno_list) <- c("bam_file")

meta <- read.csv("~/Andropogon/pca/highcov/geno_meta_01122024.csv", header = T)
str(meta)

## Add East and West tags 
east_pops <- c("TWE", "SUT", "FUL", "DES", "HAL", "MAN", "ESL", "CUI", "MON", "WEK", "WCG", "WIC", "KEN", "WAL" )
west_pops <- c("BOU","KON","MIL","AUS", "SAL","CDB","WEB","REL","BAR","AFT" )

meta$ggroup <- lapply(meta$pop_short, function(x) ifelse(x %in% east_pops, "East", "West")) %>% unlist

# Make data frame ----
df <- merge(meta, geno_list, "bam_file") %>% cbind(extracted_data_list)
str(df)
head(df)
dim(df)

## Drop clones ----
# meta[meta$pop_short == 'AFT',]
df_noclones <- df[!df$geno_short == 'Afton_TX_1',]
dim(df_noclones)

## Make long
df_long <- gather(df_noclones, chr, fis, paste(chr)) # long dataframe
df_long[,1:5] <- lapply(df_long[,1:5], factor)
str(df_long)

# Mean Fis ----
avg_inbreeding <- mean(rowMeans(df_noclones[,6:35]))
fis <- cbind(rowMeans(df_noclones[,6:35]), rowSds(as.matrix(df_noclones[,6:35])))
colnames(fis) <- c("mean", "sd")
rownames(fis) <- df_noclones$geno_short 

# plot
as.data.frame(fis) %>%
  ggplot(aes(x = reorder(rownames(fis), mean), y = mean)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Export
# write.csv(fis, file = paste0("~/Andropogon/inbreeding/mean_fis_values.",Sys.Date(),".csv"), row.names = T, quote = F)

# by chromosome
summary(df_noclones[,6:35]) 

# Plot by pop ----
# group_by(df_long, pop_short) %>%
ggplot(df_long, aes(x = geno_short, y = fis, fill = pop_short)) +
  geom_hline(yintercept = avg_inbreeding, color = "gray", linetype="dashed") +
  geom_violin() +
  # geom_jitter(size = 1) + 
  # facet_wrap(.~ggroup, scales = "free", nrow = 2) +
  # facet_wrap(.~pop_short, scales = "free", nrow = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        # legend.position = "none"
        ) +
  scale_y_continuous(limits = c(0,0.75)) +
  xlab("Genotype") + ylab("Inbreeding coefficient") +
  labs(fill="Population")

ggsave(paste0("~/Andropogon/inbreeding/inbreeding_per_genotype.", Sys.Date(),".jpeg" ),
       dpi = 300,
       height = 5,
       width = 11,
       unit = 'in')

# Plot by chrom ----
ggplot(df_long, aes(x = chr, y = fis)) +
  geom_boxplot() +
  # geom_jitter(size = 1) + 
  # facet_wrap(.~chr, scales = "free", nrow = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Model ----
mod <- lm(fis ~ pop_short + geno_short, df_long)
inbred_means <- emmeans(mod, specs = 'pop_short')
cld(inbred_means, level = 0.05)
