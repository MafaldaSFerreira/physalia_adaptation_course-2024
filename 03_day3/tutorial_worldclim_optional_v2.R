# (Optional) Script to extract abiotic environmental variables for each location 
# of interest from WorldClim
# Initial version developed by Claire et al. 2023
# Updated on June 28, 2024


# Prepare the working environment -----------------------------------------

# If not done already, install the geodata package (by executing the line below)
#install.packages("geodata")

# Load packages
library(raster)
library(geodata)

# set working directory
setwd("~/Desktop/03_day3")

# Load data ---------------------------------------------------------------

# load location data
location_GPS <- read.delim("02_data/info_pop_geo_eco.txt")

# load environmental data from WorldClim

# Ee will use the worldclim_global() function from the geodata package to 
# download the bioclimatic variable data from WorldClim.

# This function requires three these information:
#  var = "bio": This tells worldclim_global() that we want to download all 19 of the bioclimatic variables, rather than individual temperature or precipitation measurements. 
#2. res = 2.5: This is the resolution of the data we want to download; in this case, it is 2.5
# minutes of a degree. For other resolutions, you can check the documentation by typing
# ?worldclim_global into the console.
# 3. path = "data/": Finally, this sets the location to which the files are downloaded. In our case, it is the data folder we created at the beginning.

# Note also that after the files are downloaded to the data folder, they are read into memory
# and stored in the variable called bioclim_data.
r <- worldclim_global(var = "bio", res = 2.5, path = "02_data/") # this could take a few minutes
div = 10 # precision of the data

# 1 is mean temp, 12 is annual precipitations, etc...
Annual_mean_temp <- r[[1]]
variable <- paste0("bio1")


# Make a plot of the area -------------------------------------------------

# get the geographic extent area from the location data of our samples
geographic_extent_area <- extent(min(location_GPS$longitude) - 1, 
                                 max(location_GPS$longitude) + 0.5, 
                                 min(location_GPS$latitude) - 1, 
                                 max(location_GPS$latitude) + 0.5)
# plot the base map
plot((crop(Annual_mean_temp, geographic_extent_area) / div))
# add the point for the locations
points(location_GPS$longitude, location_GPS$latitude, pch = 19, col = 1, cex = 2)

### Example of how to get env data around a location point of your choice (e.g. pop L in out locations data file)
i = 1
# determine the coordinates around your point
long_min <- floor(location_GPS$longitude[i] * 10) / 10
long_max <- ceiling(location_GPS$longitude[i] * 10) / 10
lat_min <- floor(location_GPS$latitude[i] * 10) / 10
lat_max <- ceiling(location_GPS$latitude[i] * 10) / 10

# prepare the area
geographic_extent <- extent(long_min, long_max, lat_min, lat_max)

# get the value of the layer in the area
Annual_mean_temp.crop <- crop(Annual_mean_temp, geographic_extent)
mean_value_i <- global(Annual_mean_temp.crop, mean, na.rm = TRUE) / div
#mean_value_i <- mean(Annual_mean_temp.crop, na.rm = TRUE) / div
range_value_i <- (range(Annual_mean_temp.crop, na.rm = TRUE)[2] - range(Annual_mean_temp.crop, na.rm = TRUE)[1]) / div
#range_value_i <- (range(Annual_mean_temp.crop, na.rm = TRUE)[2] - range(Annual_mean_temp.crop, na.rm = TRUE)[1]) / div

# print value
location_GPS[i, ]
mean_value_i
range_value_i

##### How to plot all bioclim variables from 1 to 12
par(mfrow = c(4, 3))

for(j in 1:12) {
  x <- r[[j]] 
  div = 10
  variable <- c(variable, paste0("bioclim", j))
  geographic_extent_area <- extent(min(location_GPS$longitude) - 1, 
                                   max(location_GPS$longitude) + 0.5, 
                                   min(location_GPS$latitude) - 1, 
                                   max(location_GPS$latitude) + 0.5)

  plot((crop(x,geographic_extent_area) / div), main = variable[j])
  points(location_GPS$longitude, location_GPS$latitude, pch = 19, col = 1, cex = 2)
  
}

#### How to record all bioclim var from 1 to 12

# initialize matrix
CLIM <- matrix(nrow = dim(location_GPS), ncol = 12)
CLIM_range <- matrix(nrow = dim(location_GPS), ncol = 12)
variable <- vector(length = 0)

# loop over variables
for(j in 1:12) {
  #j = 1
  
  # get layer of variable j
  x <- r[[j]] 
  div = 10
  # name it
  variable <- c(variable, paste0("bioclim", j))
  # plot it
  geographic_extent_area <- extent(min (location_GPS$longitude) - 1, 
                                   max(location_GPS$longitude) + 0.5, 
                                   min(location_GPS$latitude) - 1, 
                                   max(location_GPS$latitude) + 0.5)

  plot((crop(x, geographic_extent_area) / div), main = variable[j])
  points(location_GPS$longitude, location_GPS$latitude, pch = 19, col = 1, cex = 2)
    
  # initialize vector for locations
  clim <- vector(length = 0)
  clim_range <- vector(length = 0)
  
  # loop over populations
  for(i in 1:dim(location_GPS)[1]) {
    #i = 1
    
    # get square around pop i
    long_min <- floor(location_GPS$longitude[i] * 10) / 10
    long_max <- ceiling(location_GPS$longitude[i] * 10) / 10
    lat_min <- floor(location_GPS$latitude[i] * 10) / 10
    lat_max <- ceiling(location_GPS$latitude[i] * 10) / 10
    geographic_extent <- extent(long_min, long_max, lat_min, lat_max)
  
    # extract layer around pop i
    x.crop <- crop(x, geographic_extent)
  
    # store info for variable j in a vector
    clim <- c(clim, global(x.crop, mean, na.rm = TRUE) / div)
    #clim <- c(clim, mean(x.crop, na.rm = TRUE) / div)
    clim_range <- c(clim_range, (range(x.crop, na.rm = TRUE)[2] - range(x.crop, na.rm = TRUE)[1]) / div)
    #clim_range <-c (clim_range, (range(x.crop, na.rm = TRUE)[2] - range(x.crop, na.rm = TRUE)[1]) / div)

  }#End for loop over populations
  
  # store info for variable j and all pop in a matrix
  #j = 1
  CLIM[, j] <- clim
  CLIM_range[, j] <- clim_range
}

colnames(CLIM) <- variable
colnames(CLIM_range) <- variable
location_CLIM <- cbind(location_GPS, CLIM)
location_CLIM_range <- cbind(location_GPS, CLIM_range)
print(location_CLIM_range)

write.table(location_CLIM, file = "02_data/worldclim.txt", quote = FALSE)
print (location_CLIM_range)



#### OCEAN CLIMATE LAYERS FOR MARINE SPATIAL ECOLOGY -- TOOLS FOR NICHE MODELING	####

#### esrigrid2ascii function - written by Elizabeth J. Sbrocco, National Evolutionary	####
#### Synthesis Center, Durham, NC  27705 elizabeth.sbrocco@gmail.com - August 2012	####

#### THE FOLLOWING IS A CUSTOM FUNCTION THAT TAKES ADVANTAGE OF THE RASTER PACKAGE TO	####
#### CROP MARSPEC ESRI GRIDS TO YOUR AREA OF INTEREST AND CONVERT IT TO AN ASCII FILE	####
#### FORMAT COMPATIBLE WITH THE MAXENT ECOLOGICAL NICHE MODELING PROGRAM.  AT THIS TIME	####
#### I HAVE ONLY WRITTEN CODE TO HANDLE ONE RASTER AT A TIME, BUT I AM WORKING ON CODE	####
#### TO DO BATCH CONVERSIONS.   NOTE THIS FUNCTION DEPENDS ON THE RASTER	####
#### PACKAGE, WHICH IN TURN DEPENDS ON THE SP PACKAGE, SO BE SURE THAT THESE ARE 	####
#### DOWNLOADED BEFORE BEGINNING.							####

#### INSTRUCTIONS:
#### AFTER DOWNLOADING THE MARSPEC DATA, UNZIP THE FILES USING 7-ZIP (OR YOUR FAVORITE	####
#### ZIP ARCHIVE UTILITY)AND MOVE THE FILES TO YOUR WORKING DIRECTORY. OPEN R AND USE 	####	
#### "FILE->CHANGE DIR" TO SELECT THE DIRECTORY CONTAINING YOUR MARSPEC GRID FILES.	####
#### EXECUTE THE FOLLOWING CODE BY COPYING AND PASTING INTO THE R CONSOLE WINDOW, OR 	####
#### HIGHLIGHT AND PRESS CTRL-R								####

##Five bioclimatic MARSPEC layers were calculated from the scaled 
#climatologies: mean annual sea surface salinity (biogeo08), 
#salinity of the least salty month (biogeo09), 
#salinity of the saltiest month (biogeo10), 
#annual range in sea surface salinity (biogeo11), 
#and annual variance in sea surface salinity (biogeo12).
#mean annual sea surface temperature (biogeo13), 
#temperature of the coldest ice-free month (biogeo14), 
#temperature of the warmest ice-free month (biogeo15), 
#annual range in sea surface temperature (biogeo16), 
#and annual variance in sea surface temperature (biogeo17).
####

##Using data from MERSPEC, we raster and select the value within a 0.1?? of latitude square around 
##the location of interest. Then within that square, I take the mean (as the value I want to use)
## but alos look at the range within the square. If it is very large, it is probably not safe 
##to do the mean on the square
##some location will find no data within that square

setwd("MARSPEC/biogeo08_17_30s") # go into the dowloaded folder with MARSPEC data to the 30s precision


# run nearly the same  code as worldlim
d <- c("08","09","10","11","12","13","14","15","16","17") # 08=mean SSS sea surface salinity, 11 = range SSS
#12=sd SSS, 13=mean SST (temp), 16=range SSt, 17=sd SST


Mean_SSS <- raster(paste("biogeo08_30s", sep = ""))
variable <- paste("biogeo08", sep = "")

# make a plot of the area
geographic_extent_area <- extent(min(location_GPS$longitude) - 1, 
max(location_GPS$longitude) + 0.5, 
min(location_GPS$latitude) - 1, 
max(location_GPS$latitude) + 0.5)

plot((crop(Mean_SSS, geographic_extent_area) / div))
points(location_GPS$longitude, location_GPS$latitude, pch = 19, col = 1, cex = 2)

# to get data round a point of your choice like pop 1
i = 1
# determine the coordinates around your point
long_min <- floor(location_GPS$longitude[i]*10)/10
long_max <- ceiling(location_GPS$longitude[i]*10)/10
lat_min <- floor(location_GPS$latitude[i]*10)/10
lat_max <- ceiling(location_GPS$latitude[i]*10)/10

# prepare the area
geographic_extent <- extent(long_min, long_max, lat_min, lat_max)

#get the value of the layer in the area
Annual_mean_temp.crop <- crop(Annual_mean_temp,geographic_extent)
mean_value_i <- mean(Annual_mean_temp.crop, na.rm = TRUE) / div
range_value_i <- (range(Annual_mean_temp.crop, na.rm = TRUE)[2] - range(Annual_mean_temp.crop, na.rm = TRUE)[1]) / div

# print value
location_GPS[i, ]
mean_value_i
range_value_i

# to plot and all  var from 1 to 10
par(mfrow = c(2,5))

for (j in 1:10) {
  x <- raster(paste("biogeo", d[j], "_30s", sep = ""))
  variable <- c(variable,paste0("marspec", d[j]))
  
  geographic_extent_area <- extent(min(location_GPS$longitude) - 1, 
  max(location_GPS$longitude) + 0.5,
  min(location_GPS$latitude) - 1, max(location_GPS$latitude) + 0.5)

  plot((crop(x, geographic_extent_area) / div), main = variable[j])
  points(location_GPS$longitude, location_GPS$latitude, pch = 19, col = 1, cex = 2)
  
}

# to record all marspec from 08 to 17
# initialize matrix:
CLIM <- matrix(nrow = dim(location_GPS), ncol = 10)
CLIM_range <- matrix(nrow = dim(location_GPS), ncol = 10)
variable <- vector(length = 0)

# loop over variables
for (j in 1:10) {
  # get layer of variable j
  x <- raster(paste("biogeo", d[j], "_30s", sep = ""))
  # name it
  variable <- c(variable, paste0("marspec",d[j]))
  # plot it
  geographic_extent_area <- extent(min(location_GPS$longitude) - 1, 
  max(location_GPS$longitude) + 0.5, 
  min(location_GPS$latitude) - 1, 
  max(location_GPS$latitude) + 0.5)

  plot((crop(x,geographic_extent_area) / div), main = variable[j])
  points(location_GPS$longitude, location_GPS$latitude, pch = 19, col = 1, cex = 2)
  
  
  # initialize vector for locations
  clim <- vector(length = 0)
  clim_range <- vector(length = 0)
  
  # loop over populations
  for (i in 1:dim(location_GPS)[1]){ 
    # get square around pop i
    long_min <- floor(location_GPS$longitude[i]*10)/10
    long_max <- ceiling(location_GPS$longitude[i]*10)/10
    lat_min <- floor(location_GPS$latitude[i]*10)/10
    lat_max <- ceiling(location_GPS$latitude[i]*10)/10
    geographic_extent <- extent(long_min, long_max, lat_min, lat_max)
    
    # extract layer around pop i
    x.crop <- crop(x, geographic_extent)
    
    # store info for variable j in a vector 
    clim <- c(clim, mean(x.crop, na.rm = TRUE) / div)
    clim_range <- c(clim_range, (range(x.crop, na.rm = TRUE)[2] - range(x.crop, na.rm = TRUE)[1]) / div)
  }
  # sotre info for variable j and all pop in a matrix
  CLIM[, j] <- clim
  CLIM_range[, j] <- clim_range
}

colnames(CLIM) <- variable
colnames(CLIM_range) <- variable
location_CLIM <- cbind(location_GPS, CLIM)
location_CLIM_range <- cbind(location_GPS, CLIM_range)
print (location_CLIM_range)

write.table(location_CLIM, file = "02_data/marspec.txt", quote = FALSE)
print(location_CLIM_range)
