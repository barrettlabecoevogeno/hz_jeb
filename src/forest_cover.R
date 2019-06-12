# Header ------------------------------------------------------------------
# Analysis of forest cover in relation to hybrid zone movement

# Load packages and functions----------------------------------------------
## Clear the environment and load in necessary packages
rm(list=ls())
library(tidyverse)
library(rethinking)
library(raster)
source("src/hzar_functions.R")


# Load data ---------------------------------------------------------

# Transects
# The empirical joint transect
joint <- read.csv("processed_data/joint_transect.csv")
# generated sites
gen <- read.csv("processed_data/generated_sites_along_transect.csv", stringsAsFactors = F)

# Rasters
# NDVI for 2000 and 2017
ndvi2000 <- raster("raw_data/forest_cover/panama_NDVI_2000_merged_clipped.tif")
ndvi2017 <- raster("raw_data/forest_cover/panama_NDVI_2017_merged_clipped.tif")
# Change in NDVI, calculated in QGIS 
# as NDVI2017-NDVI2000
dNDVI <- raster("raw_data/forest_cover/panama_deltaNDVI_2000-2017.tif")
# Lossyear for 2017
lossyear <- raster("raw_data/forest_cover/panama_lossyear_2017_merged_clipped.tif")

# Best fit cline models form 2000 and 2015
# Blum
load("results/blum_bestCline.Rdata")
blum <- right
rm(right)

# Thurman
load("results/thurman_bestCline.Rdata")
thurman <- right
rm(right)


# Process data ------------------------------------------------------------
# Make SpatialPoints objects of the collecting sites
# and of generated points along the transect.
# These spatial points objects use the same
# CRS as the Hansen data (taken from the dNDVI)
sites_actual <- SpatialPoints(cbind(joint$coord.W.decdeg, joint$coord.N.decdeg), 
                                 proj4string = CRS(proj4string(dNDVI)))

sites_gen<- SpatialPoints(cbind(gen$tran.Coord.W, gen$tran.Coord.N), 
                                 proj4string = CRS(proj4string(dNDVI)))

# Extract average NDVI and dNDVI from the rasters.
# Tested various radii, all were highly correlated.
# Chose 5k in the end, a scale that works well given our 
# distances between generated sites. 
gen$ndvi2000 <- raster::extract(ndvi2000, sites_gen, method = "simple", fun = mean, buffer = 5000, na.rm = T)
gen$ndvi2017  <- raster::extract(ndvi2017, sites_gen, method = "simple", fun = mean, buffer = 5000, na.rm = T)
gen$dNDVI<- raster::extract(dNDVI, sites_gen, method = "simple", fun = mean, buffer = 5000, na.rm = T)

# From the Lossyear file, can calculate the proportion of
# forest cover lost.
# The data is coded as 0 (no loss), and 1-17 (the year of primary loss).
# Thus, could take the proportion of forest loss (1- proportion of 0s)
prop_lost <- function(vector, ...) {
  1- (sum(vector == 0)/length(vector))
}
gen$propLost<- raster::extract(lossyear, sites_gen, method = "simple", fun = prop_lost, buffer = 5000, na.rm = T)


# Get predicted allele frequencies for generated sites from the clines
# Get the values for the cline equations
precis(blum)
precis(thurman)

# For now, will just get the point estimates, but could
# go full Bayesian and get the distribution of the allele frequency differences. 
gen$exp.p.blum <- right.eqn.vec(transectDist = gen$transect.dist, 
                                center = 467.15, width = 59.52, 
                                pmin = 0.04, pmax = 0.95, 
                                deltaR = 20.56, tauR = 0.57)
gen$exp.p.thurm <- right.eqn.vec(transectDist = gen$transect.dist, 
                                 center = 450.89, width = 93.18, 
                                 pmin = 0.04, pmax = 0.94, 
                                 deltaR = 24.81, tauR = 0.66)
gen$exp.delta.p <- gen$exp.p.thurm - gen$exp.p.blum 




# Analysis ----------------------------------------------------------------

# Correlation between
# the proportion of forest loss
# and the expected change in allele frequency
cor.test(gen$propLost, gen$exp.delta.p, method = "spearman")


# Correlation between change in NDVI and change in allele frequency
cor.test(gen$dNDVI, gen$exp.delta.p, method = "spearman")
# There is none.
