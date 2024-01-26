### Calculating Growing Degree Days using Daymet
### Author: Britney Solomon
### Date: 5/1/23

# install.packages("daymetr")
# install.packages("Rcpp")

library(daymetr) # V 1.7
library(Rcpp) # V 1.0.10
library(dplyr) # V 1.0.10

### Input climate data ----
## Common garden
# latlong <- read.csv('~/Documents/R library/2023 Stomata Project/Raw Data/andro_sitelatlong.csv', header = F)
# colnames(latlong) <- c("pop", "lat", "long")
# latlong$pop[c(2)] <- "BOU2" # Rename BOU population to prevent files from overwriting in for loop 

## All Andropogon
latlong <- read.csv("~/Andropogon/bedassle/pop_lat_long.csv", header = T)
colnames(latlong) <- c("pop", "pop2", "lat", "long", 'el')
dim(latlong)

# Create empty list
all_gdd <- list()

# Calculate GDD using a for loop ---- 
dir = "~/Andropogon/bedassle/Daymet_GDD/"

for (i in 1:length(latlong$pop)){
  
  # Create directory for each population
  dpath <- paste0(dir, latlong$pop[i], "/")
  dir.create(dpath)
  
  # Use Daymet to obtain min and max daily temperature from year range 1980-2000
  download_daymet(site = latlong$pop[i], 
                       lat = latlong$lat[i],
                       lon = latlong$long[i],
                       start = 1980, end = 2000, 
                       path = dpath, internal = FALSE)
  
  # Read Daymet data
  clm <- read.csv(paste0(dpath, latlong$pop[i], "_1980_2000.csv"), header = T, skip = 7)
  colnames(clm)[7:8] <- c("tmax", "tmin") # Simplify column names
 
  # Calculate daily GDD
  clm$tavg <- ((clm$tmax + clm$tmin)/2) # Calculate average daily temperature
  clm$gdd <- (clm$tavg - 10) # Subtract by base temp (10C)
  clm$gdd[clm$gdd < 0] <- 0 # Replace negative values with 0
  
  # Save as a file
  fileclm <- paste0(dpath, latlong$pop[i], "_GDD.csv")
  write.csv(clm, file = fileclm, row.names = F)

  # Sum GDD by year
  gdd_year = aggregate(x = clm$gdd, list(clm$year), FUN = sum)
  
  
  # Fill empty list with calculated GDD for each population
  all_gdd[[i]] <- gdd_year$x
  
}

# Turn list into data frame -- data frame contains yearly GDD for each population
yearly_gdd <- as.data.frame(do.call(cbind, all_gdd))
colnames(yearly_gdd) <- latlong$pop # Rename columns by population code
yearly_gdd$year <- gdd_year$Group.1 # Add year column
yearly_gdd <- yearly_gdd[, c(length(latlong$pop) + 1, 1:length(latlong$pop))] # Reorder dataframe so year column is first
head(yearly_gdd)

# write.csv(yearly_gdd, file = "~/Documents/R library/2023 Stomata Project/Raw Data/andro_yearly_gdd.csv" , row.names = F)


### Mean yearly GDD by population ----
avg_gdd <- colMeans(yearly_gdd)[-1]
avg_gdd
length(avg_gdd)

write.csv(avg_gdd, file = paste0(dir, "allandro_avg_gdd.csv") , row.names = T, quote = F)


