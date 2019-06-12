############################################
## Tim Thurman                            ##
## Code for analysis of                   ##
## H. erato hybrid zone in eastern Panama ##
## Constructing figure 1                  ##
############################################

# This just makes the map for figure 1,
# the rest of the figure was assembled in inkscape


# Load packages and functions----------------------------------------------
## Clear the environment and load in necessary packages
rm(list=ls())
library(tidyverse)
library("maps")
library("mapplots")
library("mapdata")
library("GISTools")
library(ggmap)
library(rgdal)
library(raster)
library(broom)
source("src/functions.R")



# Load data ---------------------------------------------------------------
joint <- read.csv("processed_data/joint_transect.csv")
ours <- joint %>%
  filter(year == "2015")
cubic.model <- lm(coord.N.decdeg ~ coord.W.decdeg + I(coord.W.decdeg^2) + I(coord.W.decdeg^3), data = joint)
cm.table <- tidy(lm(coord.N.decdeg ~ coord.W.decdeg + I(coord.W.decdeg^2) + I(coord.W.decdeg^3), data = joint))


colors <- c("darkorange", "orchid3", "green", "dodgerblue3")
# order for color and shape is blum, mallet, us
col2 <- c(alpha("grey80", 0.5), alpha("grey80", 0.5), "gold")
shapes <- c(22, 24, 21)
x.pos <-  c(-81.50000, -81.20625, -80.91250, -80.61875, -80.32500, -80.03125, -79.73750, -79.44375, -79.15000, -78.85625, -78.56250, -78.26875, -77.97500, -77.68125, -77.38750, -77.2, -77.09375)
y.pos.2 <- c(9.205280, 9.354068, 9.476545, 9.572710, 9.642563, 9.686105, 9.703334, 9.694253, 9.658859, 9.597154, 
             9.509137, 9.394809, 9.254168, 9.087216, 8.893953, 8.65, 8.45)
par(cex = 1)

#pdf("results/figure_1_map.pdf", width= 30, height=20)
par(cex = 4)
x <- maps::map("worldHires", "Panama", 
               xlim = c(-83, -76.75), ylim = c(7.2,10), 
               fill = T, col = alpha("#0e4711", alpha = 0.5), lwd = 4)
water <- readOGR("raw_data/canal_shapes/PAN_water_areas_dcw.shp")
pie.x <- seq(from = -81.5, -76.8, length.out = dim(ours)[1])
plot(water, add =T, col = "lightblue", fill = T)
curve(cm.table$estimate[1] + cm.table$estimate[2]*x + 
        cm.table$estimate[3]*x^2 ++ cm.table$estimate[4]*x^3, 
      from = -82.9, to = -77.5, add = T, lwd = 5, lty = 2)
for (row in 1:dim(ours)[1]) {
  segments(x0 = ours$coord.W.decdeg[row],
           x = x.pos[row], 
           y0 = ours$coord.N.decdeg[row],
           y = y.pos.2[row] + 0.1, lwd= 4, col = "grey40", lty = 3)
  add.pie(z = c(ours$A.melanized[row], ours$B.hetero[row], 
                ours$C.west.col[row], ours$D.postman[row]), 
          x = x.pos[row], y = y.pos.2[row] + 0.1,
          radius = 0.085, labels = NA, col = colors)
}
i <- 1
for (year1 in unique(joint$year)) {
  x <- subset(joint, year == year1)
  points(x$coord.W.decdeg, x$coord.N.decdeg, pch = shapes[i], bg = col2[i], lwd = 4)
  i <- i + 1
}
scalebar(d = 120, xy = c(-79.75,7.7), type = "bar", lonlat = T, 
         below = "kilometers", label = c(0,60,120))
north.arrow(xb = -79.325, yb = 8.3, len = 0.05, lab = "N", col = "black")
#dev.off()