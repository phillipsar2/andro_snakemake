## Title: Big Bluestem Populations Climate Data & Climate NA Testing
## Author: Britney Solomon, modified by Alyssa Phillips
## Date: 1/11/2023

library(data.table) # V 1.14.6
# install.packages('~/Downloads/ClimateNAr/', repos=NULL, type='source')
library(ClimateNAr) # V 1.1.0
library(ggplot2) # V 3.4.0

###
### Testing
###

# Climate Variables for a Single Location
test <- ClimateNAr::ClimateNA_API(ClimateBC_NA='NA', c(48.98,-115.02,200), period= 'Normal_1961_1990.nrm', MSY= 'Y')

ClimateNAr::ClimateNA_API()

head(test)

###
### Climate Data for big bluestem populations
###

# Climate Variables for Multiple Locations
input_file <- '~/cg_andropogon/data/2023 Stomata Project/Raw Data/Population Data.csv'
clm <- ClimateNAr::ClimateNA_API2(ClimateBC_NA='NA', inputFile = input_file, period = 'Normal_1961_1990.nrm', MSY= 'Y')
head(clm)
# write.csv(clm, file = "~/Documents/R library/androclimatedata.csv", row.names = F, col.names = T)

# Plotting Variables
clm <- read.csv("~/cg_andropogon/scripts/Climate Data - ClimateNA/androclimatedata_popcodes.csv")
str(clm)
dim(clm)

ggd <- read.csv("~/cg_andropogon/scripts/GDD Calculations/andro_yearly_gdd.csv")
str(ggd)
dim(ggd)

avg_ggd <- colMeans(ggd[,2:17])
length(avg_ggd)

clm$GGD <- avg_ggd
str(clm)
dim(clm)

bou <- colMeans(clm[1:2,3:31])

df <- clm[2:16,]
df[1,3:31] <- bou 

# ggplot(clm, aes(x=lat, y=MAT))+
#   geom_point()
# cor(clm$lat, clm$MAT)
# ggplot(clm, aes(x=lat, y=MAP))+
#   geom_point()
# ggplot(clm, aes(x=lat, y=PAS))+
#   geom_point()

tmp <- df[3:31]
rownames(tmp) <- df$ID2
pc <- prcomp(tmp)
biplot(pc)

plot(pc$x[,1], pc$x[,2])

cor(tmp) %>% corrplot::corrplot()
