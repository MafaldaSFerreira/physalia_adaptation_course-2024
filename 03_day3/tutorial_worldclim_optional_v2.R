# (Optional) Script to extract environmental variables for each location of interest using Worldclim

# Load packages
library(raster)
#library(rgdal)
library(terra)    # to install
library(sdmpredictors)     # to install

# load location data
location_GPS <- read.delim("02_data/info_pop_geo_eco.txt")

###
# we'll use functions of the 'sdmpredictors' package to access different online datasets
pred_datasets <- list_datasets(terrestrial = TRUE, marine = TRUE)
pred_datasets
names(pred_datasets)
pred_datasets[ , 1:4]  # these are the datasets currently available for download using the 'sdmpredictors' package; you can check their URLs for more info on their contents
pred_datasets[ , c("dataset_code", "citation")]  # remember to ALWAYS cite the actual data sources, not just the package you used for downloading!

pred_layers <- list_layers(datasets = pred_datasets)
unique(pred_layers$dataset_code)
unique(pred_layers[pred_layers$dataset_code == "WorldClim", ]$name)  # example of terrestrial variables dataset
unique(pred_layers[pred_layers$dataset_code == "MARSPEC", ]$name)  # example of marine variables dataset

# let's choose one dataset (e.g. WorldClim) and one particular set of variables (e.g. altitude and the bioclimatic ones, which are in rows 1 to 20):
layers_choice <- unique(pred_layers[pred_layers$dataset_code == "WorldClim", c("name", "layer_code")])
layers_choice
layers_choice <- layers_choice[1:20, ]
layers_choice

# define folder for downloading / fetching the variables' map layers:
options(sdmpredictors_datadir = "../outputs/sdmpredictors")
# load the layers to the current R session (downloading them if they aren't already in the folder defined above):
layers <- load_layers(layers_choice$layer_code, rasterstack = FALSE)  # rasterstack=TRUE gives error when there are layers with different extent
layers  # a list of raster maps
# see how many elements in 'layers':
length(layers)

# plot a couple of layers to see how they look:
names(layers)
# convert each layer to 'SpatRaster' class (from package 'terra'), which is much faster to process:
layers <- lapply(layers, rast)
plot(layers[[1]], main = names(layers)[1])  # ignore errors if map is plotted anyway
plot(layers[[5]], main = names(layers)[5])

# find out if your layers have different extents or resolutions:
unique(pred_layers[pred_layers$dataset_code == "WorldClim", ]$cellsize_lonlat)  # in this case 0.08333333 - spatial resolution can then be coarsened as adequate for your species data and modelling region (see below)
sapply(layers, ext)  # if you get different extents (which doesn't happen with WorldClim, but may happen with other datasets), you'll have to crop all layers to the minimum common extent before proceeding
# for example, if the first layer has the smallest extent:
#layers <- lapply(layers, crop, extent(layers[[1]]))

# once all layers have the same extent and resolution, you can combine them in a single multi-layer raster map:
layers <- rast(layers)
layers
plot(layers)







###

r <- raster::getData("worldclim", var = "bio", res = 2.5)
div = 10 # precision of the data

# 1 is mean temp, 12 is annual precipitations, etc...
Annual_mean_temp <- r[[1]]
variable <- paste0("bio1")

# make a plot of the area
aoi_area <- extent(min(location_GPS$longitude) - 1, 
max(location_GPS$longitude) + 0.5, 
min(location_GPS$latitude) - 1, 
max(location_GPS$latitude) + 0.5)

plot((crop(Annual_mean_temp, aoi_area) / div))
points(location_GPS$longitude, location_GPS$latitude, pch = 19, col = 1, cex = 2)

# to get data round a point of your choice like pop 1
i = 1
# determine the coordinates around your point
long_min <- floor(location_GPS$longitude[i] * 10) / 10
long_max <- ceiling(location_GPS$longitude[i] * 10) / 10
lat_min <- floor(location_GPS$latitude[i] * 10) / 10
lat_max <- ceiling(location_GPS$latitude[i] * 10) / 10

# prepare the area
aoi <- extent(long_min, long_max, lat_min, lat_max)

# get the value of the layer in the area
Annual_mean_temp.crop <- crop(Annual_mean_temp, aoi)
mean_value_i <- mean(Annual_mean_temp.crop@data@values, na.rm = TRUE) / div
range_value_i <- (range(Annual_mean_temp.crop@data@values, na.rm = TRUE)[2] - range(Annual_mean_temp.crop@data@values, na.rm = TRUE)[1]) / div

# print value
location_GPS[i, ]
mean_value_i
range_value_i

# to plot and all bioclim var from 1 to 12
par(mfrow = c(4,3))

for (j in 1:12) {
  x <- r[[j]] 
  div = 10
  variable <- c(variable, paste0("bioclim", j))
  aoi_area <- extent(min(location_GPS$longitude) - 1, 
  max(location_GPS$longitude) + 0.5, 
  min(location_GPS$latitude) - 1, 
  max(location_GPS$latitude) + 0.5)

  plot((crop(x,aoi_area) / div), main = variable[j])
  points(location_GPS$longitude, location_GPS$latitude, pch = 19, col = 1, cex = 2)
  
}

# to record all bioclim var from 1 to 12
# initialize matrix:
CLIM <- matrix(nrow = dim(location_GPS), ncol = 12)
CLIM_range <- matrix(nrow = dim(location_GPS), ncol = 12)
variable <- vector(length = 0)

# loop over variables
for (j in 1:12) {
  # get layer of variable j
  x <- r[[j]] 
  div = 10
  # name it
  variable <- c(variable, paste0("bioclim", j))
  # plot it
  aoi_area <- extent(min (location_GPS$longitude) - 1, 
  max(location_GPS$longitude) + 0.5, 
  min(location_GPS$latitude) - 1, 
  max(location_GPS$latitude) + 0.5)

  plot((crop(x, aoi_area) / div), main = variable[j])
  points(location_GPS$longitude, location_GPS$latitude, pch = 19, col = 1, cex = 2)
    
  # initialize vector for locations
  clim <- vector(length = 0)
  clim_range <- vector(length = 0)
  
  # loop over populations
  for (i in 1:dim(location_GPS)[1]) { 
    # get square around pop i
    long_min <- floor(location_GPS$longitude[i] * 10) / 10
    long_max <- ceiling(location_GPS$longitude[i] * 10) / 10
    lat_min <- floor(location_GPS$latitude[i] * 10) / 10
    lat_max <- ceiling(location_GPS$latitude[i] * 10) / 10
    aoi <- extent(long_min, long_max, lat_min, lat_max)
  
    # extract layer around pop i
    x.crop <- crop(x, aoi)
  
    # store info for variable j in a vector 
    clim <-c (clim, mean(x.crop@data@values, na.rm = TRUE) / div)
    clim_range <-c (clim_range, (range(x.crop@data@values, na.rm = TRUE)[2] - range(x.crop@data@values, na.rm = TRUE)[1]) / div)

  }
  # store info for variable j and all pop in a matrix
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
aoi_area <- extent(min(location_GPS$longitude) - 1, 
max(location_GPS$longitude) + 0.5, 
min(location_GPS$latitude) - 1, 
max(location_GPS$latitude) + 0.5)

plot((crop(Mean_SSS, aoi_area) / div))
points(location_GPS$longitude, location_GPS$latitude, pch = 19, col = 1, cex = 2)

# to get data round a point of your choice like pop 1
i = 1
# determine the coordinates around your point
long_min <- floor(location_GPS$longitude[i]*10)/10
long_max <- ceiling(location_GPS$longitude[i]*10)/10
lat_min <- floor(location_GPS$latitude[i]*10)/10
lat_max <- ceiling(location_GPS$latitude[i]*10)/10

# prepare the area
aoi <- extent(long_min, long_max, lat_min, lat_max)

#get the value of the layer in the area
Annual_mean_temp.crop <- crop(Annual_mean_temp,aoi)
mean_value_i <- mean(Annual_mean_temp.crop@data@values, na.rm = TRUE) / div
range_value_i <- (range(Annual_mean_temp.crop@data@values, na.rm = TRUE)[2] - range(Annual_mean_temp.crop@data@values, na.rm = TRUE)[1]) / div

# print value
location_GPS[i, ]
mean_value_i
range_value_i

# to plot and all  var from 1 to 10
par(mfrow = c(2,5))

for (j in 1:10) {
  x <- raster(paste("biogeo", d[j], "_30s", sep = ""))
  variable <- c(variable,paste0("marspec", d[j]))
  
  aoi_area <- extent(min(location_GPS$longitude) - 1, 
  max(location_GPS$longitude) + 0.5,
  min(location_GPS$latitude) - 1, max(location_GPS$latitude) + 0.5)

  plot((crop(x, aoi_area) / div), main = variable[j])
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
  aoi_area <- extent(min(location_GPS$longitude) - 1, 
  max(location_GPS$longitude) + 0.5, 
  min(location_GPS$latitude) - 1, 
  max(location_GPS$latitude) + 0.5)

  plot((crop(x,aoi_area) / div), main = variable[j])
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
    aoi <- extent(long_min, long_max, lat_min, lat_max)
    
    # extract layer around pop i
    x.crop <- crop(x, aoi)
    
    # store info for variable j in a vector 
    clim <- c(clim, mean(x.crop@data@values, na.rm = TRUE) / div)
    clim_range <- c(clim_range, (range(x.crop@data@values, na.rm = TRUE)[2] - range(x.crop@data@values, na.rm = TRUE)[1]) / div)
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
