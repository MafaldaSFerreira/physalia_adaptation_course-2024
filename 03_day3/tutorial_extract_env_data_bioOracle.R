# Extract oceanographic data layers from Bio-Oracle v.2.1
# https://www.bio-oracle.org/code.php
# Angela Fuentes Pardo, e-mail: apfuentesp@gmail.com

# NOTE: This code 
# The data available in Bio-ORACLE v.2.1 are documented in two peer reviewed articles that you should cite:
  
# [1] Tyberghein L, Verbruggen H, Pauly K, Troupin C, Mineur F, De Clerck O (2012) Bio-ORACLE: A global environmental dataset 
# for marine species distribution modelling. Global Ecology and Biogeography, 21, 272–281.
# [2] Assis, J., Tyberghein, L., Bosh, S., Verbruggen, H., Serrão, E. A., & De Clerck, O. (2017). Bio-ORACLE v2.0: Extending 
# marine data layers for bioclimatic modelling. Global Ecology and Biogeography.


# Setup working environment -----------------------------------------------

# Clean environment space
rm(list=ls())

# Set environment variables
working_dir <- "~/Desktop/03_day3"
metadata_file <- "02_data/herring_locations.csv"

# Set working directory
setwd(working_dir)

# Load packages 
library(sdmpredictors)
library(leaflet)

# explore all available layers 
list_layers()

# Extract environmental information for locations -------------------------

# Download environmental data layers (here we chose 6)
env_list <- c("BO21_salinitymean_bdmean", "BO21_salinityrange_bdmean", "BO21_tempmean_bdmean", "BO21_temprange_bdmean", "BO21_dissoxmean_bdmean", "BO21_dissoxrange_bdmean")
environment.bottom <- load_layers(layercodes = c(env_list), equalarea = FALSE, rasterstack = TRUE, datadir = "02_data") 

# Download bathymetry 
bathymetry <- load_layers("BO_bathymean", datadir = "02_data") 

# Load data frame with the location information (geographic coordinates, latitude and longitude) 
my_sites <- read.table(metadata_file, header = TRUE, sep = ";", comment.char = "", stringsAsFactors = FALSE)
my_sites
str(my_sites)
dim(my_sites)

# Visualise sites of interest in google maps 
m <- leaflet(data = my_sites) %>% 
  addTiles() %>%
  addMarkers(~longitude, ~latitude, popup = ~as.character(sample), label = ~as.character(sample))
m

# Extract environmental values from layers 
env <- data.frame(Sample = my_sites$sample,
                  #Code=my_sites$Code,
                  depth = raster::extract(bathymetry, my_sites[, c("longitude", "latitude")]), 
                  raster::extract(environment.bottom, my_sites[, c("longitude", "latitude")], buffer = 10000, fun = mean) # radius of m around center point, and calculate average of all cells found
                  ) 
head(env)

# Save env data in a file
#write.table(env, "Herring_GEA_bio-oracle_SST_min_max_min_depth_62pools.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


# Examine correlation between environmental variables ---------------------

# Load packages
library(psych)    # to investigate correlations among predictors
library(data.table)

# Explore the structure of the dataframe. Make sure enviro variables are numeric (num)
str(env)  
head(env)

# Examine whether there is high correlation among predictors
#pdf('Correlation_plot_env_var.pdf', height = 15, width = 15)
pairs.panels(env[, 2:ncol(env)], scale = TRUE) # exclude first column with sample name
#dev.off()

# as you can see, there is a high correlation between (R2 > 0.7)



# Plot environmental layers with base R graphics --------------------------

# Set map boundary
geographic_extent <- raster::extent(-23, 20, 18, 69) # min Log, max Log, min Lat, max Lat

# Get list of the vriables you want to plot
env_list <- colnames(env[, 3:ncol(env)])
env_list

# Make the plots for all the variables of interest

#pdf("Maps_environmental_data_biooracle.pdf")
par(mfrow=c(3,3))
#par(mfrow=c(2,2))

for(i in env_list){
  #i="BO21_salinityrange_bdmean"
  #i="BO21_tempmean_bdmean"
  
  # Easy download of raster file 
  env_layer <- load_layers(i, datadir = "02_data")
  
  # Crop raster to fit the area of interest
  env_layer.crop <- raster::crop(env_layer, geographic_extent) 
  
  # Choose a color palette 
  library(palr)
  sst_pal <- sst_pal(palette = TRUE)
  my.colors <- colorRampPalette(sst_pal$cols) 
  
  # Plot the map 
  plot(env_layer.crop, col = my.colors(1000), axes = FALSE, box = FALSE, main = i) 
  #title(cex.sub = 1.25, sub = i) 
}

#dev.off()



# Plot environmental layers with ggplot2 ----------------------------------

remotes::install_github("ropensci/rnaturalearthhires")
install.packages("rnaturalearthhires", repos = "https://ropensci.r-universe.dev", type = "source")


# load packages
library(palr)
library(raster)
#library(sf)
library(dplyr)
library(rnaturalearth)
#library(rnaturalearthdata)
library(rnaturalearthhires)
#library(maps)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

# Set map boundary (xmin, xmax, ymin, ymax)
boundary <- extent(-23, 20, 18, 69)  # min Log, max Log, min Lat, max Lat

# Create a ggplot theme for heatmaps
ggtheme <- theme(axis.title = element_text(size = 12),
                 axis.text = element_text(size = 11, colour = "black"),
                 panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.5),
                 legend.title = element_text(size = 13),
                 legend.text = element_text(size = 11),
                 plot.title = element_text(size = 14, hjust = 0.5),
                 panel.grid = element_blank(),
                 legend.key.size = unit(1, "cm"),
                 legend.key.width = unit(0.3, 'cm'))

p_list <- list()

for(i in env_list){
  print(paste0("currently processing...", i))
  #i = "BO21_salinitymean_bdmean"

  # Easy download of raster file 
  env_layer <- load_layers(i, datadir = "02_data")
  
  # Crop raster to fit the area of interest
  env_layer.crop <- raster::crop(env_layer, boundary) %>% rasterToPoints() %>% data.frame()
  
  # Define a color palette 
  #library(palr)
  sst_pal <- sst_pal(palette = TRUE)
  my.colors <- colorRampPalette(sst_pal$cols)
  
  # Plot the map 
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
patchwork::wrap_plots(p_list, nrow = 3, ncol = 3)

# Safe the map in a file
#ggsave("Map_6_envVars_HOM.pdf", width = 10, height = 9)
#ggsave("Map_6_envVars_HOM.png", width = 10, height = 9, dpi = 300)


