# Extract environmental data for specific locations and do some data exploration
# Angela Fuentes Pardo, e-mail: apfpgen@gmail.com
# Created on: June, 2023
# Last updated on: June 28, 2024
#
# DISCLAIMER: 
# This script presents the general steps to extract environmental data available
# as raster files in a public database, such as Bio-Oracle (https://bio-oracle.org/index.php).
# Note that this code was developed for Bio-ORACLE v.2.1. On July 2023, Bio-ORACLE
# went through significant restructuring. Thus, this code will need some tweaking 
# in order to work for the latest database version v.3.0.
#
# Data in Bio-ORACLE v.2.1 are documented in two peer reviewed articles:
#  
# [1] Tyberghein L, Verbruggen H, Pauly K, Troupin C, Mineur F, De Clerck O (2012) Bio-ORACLE: A global environmental dataset 
# for marine species distribution modelling. Global Ecology and Biogeography, 21, 272–281.
# [2] Assis, J., Tyberghein, L., Bosh, S., Verbruggen, H., Serrão, E. A., & De Clerck, O. (2017). Bio-ORACLE v2.0: Extending 
# marine data layers for bioclimatic modelling. Global Ecology and Biogeography.


# (If needed) Install required R packages  --------------------------------

install.packages(c("sdmpredictors", "leaflet", "palr", "raster", "dplyr", "maps", "rnaturalearth", "ggplot2", "RColorBrewer", "patchwork"))
remotes::install_github("ropensci/rnaturalearthhires")
install.packages("rnaturalearthhires", repos = "https://ropensci.r-universe.dev", type = "source")

# Setup working environment -----------------------------------------------

# clean environment space
rm(list=ls())

# set environment variables
working_dir <- "~/Desktop/03_day3" # change as needed to match the configuration in your computer
metadata_file <- "02_data/herring_locations.csv"

# set working directory
setwd(working_dir)

# load packages 
library(sdmpredictors)
library(leaflet) # to make interactive maps
library(palr)  # color palette based on remotely sensed data
library(raster)
library(dplyr)
library(maps)
library(rnaturalearth)
library(rnaturalearthhires)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

# Download environmental data from a public database ----------------------

# explore which env variables are available in Bio-Oracle

# we will use functions of the "sdmpredictors" package to access online environmental datasets
pred_datasets <- list_datasets(terrestrial = TRUE, marine = TRUE)
pred_datasets
names(pred_datasets)
pred_datasets[, 1:4]  # these are the datasets currently available for download; you can check their URLs for more info on their contents
pred_datasets[, c("dataset_code", "citation")]  # remember to always cite the actual data sources, not just the package you used to download the data

pred_layers <- list_layers(datasets = pred_datasets)
unique(pred_layers$dataset_code)
unique(pred_layers[pred_layers$dataset_code == "WorldClim", ]$name)  # example of terrestrial variables dataset
unique(pred_layers[pred_layers$dataset_code == "MARSPEC", ]$name)  # example of marine variables dataset

# download the raster files of the variables of interest (here we chose 6)
env_list <- c("BO21_salinitymean_bdmean", "BO21_salinityrange_bdmean", "BO21_tempmean_bdmean", "BO21_temprange_bdmean", "BO21_dissoxmean_bdmean", "BO21_dissoxrange_bdmean")
env_data <- load_layers(layercodes = c(env_list), equalarea = FALSE, rasterstack = TRUE, datadir = "02_data") 

# download bathymetry data
bathymetry <- load_layers("BO_bathymean", datadir = "02_data") 

# load the location information (geographic coordinates, latitude and longitude) 
locations <- read.table(metadata_file, header = TRUE, sep = ";", comment.char = "", stringsAsFactors = FALSE)
str(locations) # explore their data types
dim(locations) # explore the dimensions of the table
head(locations)

# visualise loctions on an interactive map using leaflet 
m <- leaflet(data = locations) %>% 
  addTiles() %>%
  addMarkers(~longitude, ~latitude, popup = ~as.character(sample), label = ~as.character(sample))
m # note that with the mouse you can scrollin and out to zoom in and out in the map

# Plot environmental layers with base R graphics --------------------------

# set map boundaries
geographic_extent <- raster::extent(-70, 30, 42, 71) # min Longitude, max Longitude, min Latitude, max Latitude

# make plots

#pdf("Maps_environmental_data_biooracle.pdf")
# set plot display
par(mfrow=c(3,2))
#par(mfrow=c(3,3))

# go over the list of env variables to plot
for(i in env_list){
  #i="BO21_salinityrange_bdmean"
  #i="BO21_tempmean_bdmean"
  
  # load current raster file 
  env_layer <- load_layers(i, datadir = "02_data")
  
  # crop raster to fit the area of interest
  env_layer.crop <- raster::crop(env_layer, geographic_extent) 
  
  # choose color palette 
  sst_pal <- sst_pal(palette = TRUE)
  my.colors <- colorRampPalette(sst_pal$cols) 
  
  # plot the map 
  plot(env_layer.crop, col = my.colors(1000), axes = FALSE, box = FALSE, main = i) 
  #title(cex.sub = 1.25, sub = i) 
}#End of loop
#dev.off()


# Plot environmental layers with ggplot2 ----------------------------------

# set map boundaries
#geographic_extent <- raster::extent(-70, 30, 42, 71) # min Log, max Log, min Lat, max Lat

# set ggplot theme
ggtheme <- theme(axis.title = element_text(size = 12),
                 axis.text = element_text(size = 11, colour = "black"),
                 panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
                 legend.title = element_text(size = 13),
                 legend.text = element_text(size = 11),
                 plot.title = element_text(size = 14, hjust = 0.5),
                 panel.grid = element_blank(),
                 legend.key.size = unit(1, "cm"),
                 legend.key.width = unit(0.3, 'cm'))

# initialize an empty list where plots will be stored
p_list <- list()

# go over the list of env variables to plot
for(i in env_list) {
  print(paste0("currently processing layer : ", i))
  #i = "BO21_salinitymean_bdmean"

  # download of raster file 
  env_layer <- load_layers(i, datadir = "02_data")
  
  # crop raster to fit the area of interest
  env_layer.crop <- raster::crop(env_layer, geographic_extent) %>% rasterToPoints() %>% data.frame()
  
  # set color palette 
  sst_pal <- sst_pal(palette = TRUE)
  my.colors <- colorRampPalette(sst_pal$cols)
  
  # plot map 
  p <- ggplot() +
    geom_tile(data = env_layer.crop, aes(x = x, y = y, fill = env_layer.crop[, 3])) +
    coord_quickmap(expand = FALSE) +
    xlab("Longitude") +
    ylab("Latitude") +
    ggtitle(i) +
    scale_fill_gradientn(name = "", #expression(~degree~C), 
                         colours = my.colors(10), 
                         limits = c(min(env_layer.crop[, 3]), max(env_layer.crop[, 3])),
                         guide = guide_colorbar(frame.colour = "black", 
                                                ticks.colour = "black", linewidth = 0.1)) +
    ggtheme
  p
  #ggsave(paste0(i, ".pdf"), width = 10, height = 9)
#  ggsave(paste0(i, ".png"), width = 3, height = 4, dpi = 200)
  
  # Save plot to the list of plots
  p_list[[i]] <- p
  
}#End for loop

# combine plots
patchwork::wrap_plots(p_list, nrow = 3, ncol = 2)
#patchwork::wrap_plots(p_list, nrow = 3, ncol = 3)

# safe the map in a file
#ggsave("Map_environmental_variables.pdf", width = 10, height = 9)
#ggsave("Map_environmental_variables", width = 10, height = 9, dpi = 300)


# Extract environmental data for each location ----------------------------

# store the data into a data frame (table)
env <- data.frame(sample = locations$sample_name, #locations$sample,
                  longitude = locations$longitude,
                  latitude = locations$latitude,
                  depth = raster::extract(bathymetry, locations[, c("longitude", "latitude")]), 
                  raster::extract(env_data, locations[, c("longitude", "latitude")], 
                  # the buffer value sets the radius of m around center point, and then we calculate the mean of all cells found
                                  buffer = 5000, fun = mean) 
) 
# explore the resulting table
head(env)
# Note that for each location we have a value of each environmental variable. 
# This value corresponds to the mean calculated for a circle with radius 5 Km around the target location

# save env data in a file
write.table(env, "02_data/env_var_locations.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


# Examine an environmental variable across locations ----------------------

# let's examine how salinity varies across locations

# save in an object the name of the target env variable
target_env <- "BO21_salinitymean_bdmean"

# set color palette
myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")), bias=.5)

# load country shapes
admin <- ne_countries(scale = "medium", returnclass = "sf")

# make plot
ggplot() +
  geom_sf(data = admin, fill = "white", size = 0.1) +
  #geom_sf(data = admin, fill = gray(0.9), size = 0.1) +
  geom_point(data = env, aes(x = longitude, y = latitude, color = .data[[target_env]])) +
  scale_colour_gradientn(colours = myPalette(11), name = "salinity (PSU)", na.value = "black") +
  #scale_colour_gradient(low = "blue", high = "red", name = "salinity (PSU)", na.value = "gray") +
  guides(fill = guide_colorbar(barwidth = 0.7, 
                               barheight = 11, 
                               frame.linewidth = 0.5, 
                               ticks.linewidth = 0.5,
                               ticks.colour = "black", 
                               frame.colour="black")) +
  coord_sf(xlim = c(-70, 30), ylim = c(42, 71), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude") +
  #coord_quickmap(expand = FALSE) +
  ggtheme


# now do the same but for mean temperature
target_env <- "BO21_tempmean_bdmean"

# make plot
ggplot() +
  geom_sf(data = admin, fill = "white", size = 0.1) +
  #geom_sf(data = admin, fill = gray(0.9), size = 0.1) +
  geom_point(data = env, aes(x = longitude, y = latitude, color = .data[[target_env]])) +
  scale_colour_gradientn(colours = my.colors(1000), name = "SST (C°)", na.value = "black") +
  #scale_colour_gradientn(colours = myPalette(11), name = "SST (C°)", na.value = "black") +
  #scale_colour_gradient(low = "blue", high = "red", name = "SST (C°)", na.value = "gray") +
  #scale_colour_gradient(colours = myPalette(11), name = "SST (C°)", na.value = "gray") +
  guides(fill = guide_colorbar(barwidth = 0.7, 
                               barheight = 11, 
                               frame.linewidth = 0.5, 
                               ticks.linewidth = 0.5,
                               ticks.colour = "black", 
                               frame.colour="black")) +
  coord_sf(xlim = c(-70, 30), ylim = c(42, 71), expand = FALSE) +
  xlab("Longitude") + ylab("Latitude") +
  #coord_quickmap(expand = FALSE) +
  ggtheme


# Examine the extent of correlation between environmental variables -------

# Load packages
library(psych) # to investigate correlations among predictors
library(data.table)

# explore the structure of the dataframe. Make sure enviro variables are numeric (num)
str(env)  
head(env)

# estimate the correlation coefficient between pairs of environmental variables (predictors for RDA)
#pdf('Correlation_plot_env_var.pdf', height = 15, width = 15)
pairs.panels(env[, 2:ncol(env)], scale = TRUE) # exclude first column with sample name
#dev.off()

# as you may notice, there is a high correlation between some variables (|R2 > 0.7|),
# e.g., between temperature range and O2 dissolved range (R2 = 0.94).
# For RDA, we may need to account for this, by either, excluding one of the two 
# highly correlated variables, or by creating synthetic variables using PCA.

