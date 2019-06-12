############################################
## Tim Thurman                            ##
## Code for analysis of                   ##
## H. erato hybrid zone in eastern Panama ##
## Figure 3                               ##
############################################

# Load packages -----------------------------------------------------------
# load packages
library(tidyverse)
library(cowplot)
library(bayesplot)
library(raster)
library(broom)
source("src/hzar_functions.R")


# library(rlang)
# 
# library(assertthat)
# library(coda)
# library(rstan)
# library(kableExtra)
# library(rethinking)
# library(rgdal)
# 
# library(rgeos)
# library(sp)
# library(sf)
# 
# library(stringr)
# library(ggrepel)

# source("src/functions.R")


# Load data ---------------------------------------------------------------
# Read in the generated sites along the cline
gen <- read.csv("processed_data/generated_sites_along_transect.csv", stringsAsFactors = F)
# read in the data on change in NDVI
dNDVI <- raster("raw_data/forest_cover/panama_deltaNDVI_2000-2017.tif")
# Turn sites into spatial points objects
sites_gen<- SpatialPoints(cbind(gen$tran.Coord.W, gen$tran.Coord.N), 
                          proj4string = CRS(proj4string(dNDVI)))
# and calculate the mean dNDVI around each generated site, 5k radius
dNDVI_gen_5k <- raster::extract(dNDVI, sites_gen, method = "simple", fun = mean, buffer = 5000, na.rm = T)
gen$dNDVI <- dNDVI_gen_5k

# Lossyear for 2017
lossyear <- raster("raw_data/forest_cover/panama_lossyear_2017_merged_clipped.tif")

prop_lost <- function(vector, ...) {
  1- (sum(vector == 0)/length(vector))
}
gen$propLost<- raster::extract(lossyear, sites_gen, method = "simple", fun = prop_lost, buffer = 5000, na.rm = T)

#Then, get expected change in allele frequency for each site
gen$exp.p.blum <- right.eqn.vec(transectDist = gen$transect.dist, 
                                center = 467.15, width = 59.52, 
                                pmin = 0.04, pmax = 0.95, 
                                deltaR = 20.56, tauR = 0.57)
gen$exp.p.thurm <- right.eqn.vec(transectDist = gen$transect.dist, 
                                 center = 450.89, width = 93.18, 
                                 pmin = 0.04, pmax = 0.94, 
                                 deltaR = 24.81, tauR = 0.66)
gen$exp.delta.p <- gen$exp.p.thurm - gen$exp.p.blum


r <- round(cor(gen$propLost, gen$exp.delta.p, method = "spearman"), digits = 2)
p <- round(tidy(cor.test(gen$propLost, gen$exp.delta.p, method = "spearman"))$p.value, digits = 2)

r2 <- round(cor(gen$dNDVI, gen$exp.delta.p, method = "spearman"), digits = 2)
p2 <- round(tidy(cor.test(gen$dNDVI, gen$exp.delta.p, method = "spearman"))$p.value, digits = 2)




# Plot --------------------------------------------------------------------
twoA <- ggplot(data = gen, aes(x = propLost, y = exp.delta.p)) +
  geom_point(size = 3, color = alpha("darkgreen", alpha = 0.7)) +
  draw_label(x = 0.23, y = 0.19, label = expression(rho*" = -0.31, p = 0.04"), size = 12) +
  theme_default() + 
  xlab("Proportion of forest lost, 2000-2017") +
  ylab(expression(Delta*"Cr"[HYD]*" 1999-2015")) +
  theme(axis.text = element_text(size = 10))


twoB <- ggplot(data = gen, aes(x = dNDVI, y = exp.delta.p)) +
  geom_point(size = 3, color = alpha("darkgreen", alpha = 0.7))  +
  theme_default() + 
  draw_label(x = -0.075, y = .19, label = expression(rho*" = -0.16, p = 0.27"), size = 12) +
  xlab(expression(Delta*"NDVI, 2000-2017")) +
  ylab(expression(Delta*"Cr"[HYD]*" 1999-2015")) +
  theme(axis.text = element_text(size = 10))

a <- plot_grid(twoA, twoB, labels = c("A", "B"), align = "hv", rel_widths = c(1,1), rel_heights = c(1,1), axis = "tblr")

ggsave("results/figure_3.pdf", a,  height = 83, width = 166, units = "mm", dpi = 1000)
