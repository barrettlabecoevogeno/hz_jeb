############################################
## Tim Thurman                            ##
## Code for analysis of                   ##
## H. erato hybrid zone in eastern Panama ##
## Combining three collecting periods     ##
## onto the same transect                 ##
############################################


# Load in packages --------------------------------------------------------
rm(list=ls())
library(assertthat)
library(readxl)
library("geosphere")
library("maps")
library("GISTools")
library("mapdata")
library(tidyverse)
library(broom)
source("src/functions.R")


# Load in mallet data, convert GPS_-----------------------------------------
# Mallet sampled 20 sites total, but didn't include all of them for his transect calculations.
# We won't use all sites either.
# Exclude: site 3 (on an island off Panama)
#          Sites 15 and 19: On northern coast
#          Sites 16 and 18 (It seems Mallet didn't include these either, but Blum included 18)
# We will include two sites that Mallet excluded. Sites 12 and 20. Both are quite close to the main
# transect, and can easily be included given our cubic regression framework.

# This leaves us with 15 Mallet sites to estimate the cline. 

m.main <- c(rep("Y", times = 20))
m.main[c(3,15,16,18,19)] <- "N"
mallet <- read_excel("raw_data/mallet_1986/mallet_data_table4_erato.xlsx") %>%
  separate(lat.deg.min, into = c("lat.degrees", "lat.minutes"), convert = T) %>%
  separate(long.deg.min, into = c("long.degrees", "long.minutes"), convert =T) %>%
  mutate(lat.dec = lat.minutes/60,
         long.dec = long.minutes/60) %>%
  mutate(coord.N.decdeg = lat.degrees + lat.dec,
         coord.W.decdeg = (long.degrees + long.dec)*-1) %>%
  mutate(year = 1982) %>%
  mutate(on.main = m.main) %>%
  filter(on.main == "Y") %>%
  dplyr::select(year, site.collected = site.name.full, coord.N.decdeg, coord.W.decdeg,
         A.melanized, B.hetero, C.west.col, D.postman)

# Load in blum data, convert GPS ---------------------------------------
# Blum sampled 24 sites total, but didn't include all of them for his transect calculations.
# We won't use all sites either.
# Exclude: Site 1 and Site 11 (Islands off the mainland of Panama)

# We will include two sites that Blum excluded. Sites 9 and 10 (PLR and Madden dam)
# Both are quite close to the main transect, and can easily be included given our 
# cubic regression framework. Also, we sampled at PLR, and Mallet sampled at Madden dam

# This leaves us with 22 sites for Blum.
b.main <- c(rep("Y", times = 24))
b.main[c(1, 11)] <- "N" 
blum <- read_excel("raw_data/blum_2002/blum_date_table1_erato.xlsx") %>%
  separate(lat.deg.min, into = c("lat.degrees", "lat.minutes"), convert = T, sep = " ") %>%
  separate(long.deg.min, into = c("long.degrees", "long.minutes"), convert =T, sep = " ") %>%
  mutate(lat.dec = lat.minutes/60,
         long.dec = long.minutes/60) %>%
  mutate(coord.N.decdeg = lat.degrees + lat.dec,
         coord.W.decdeg = (long.degrees + long.dec)*-1) %>%
  mutate(year = 1999) %>%
  mutate(on.main = b.main) %>%
  filter(on.main == "Y") %>%
  dplyr::select(year, site.collected = site.name.full, coord.N.decdeg, coord.W.decdeg,
         A.melanized, B.hetero, C.west.col, D.postman)
blum$C.west.col <- 0

rm(b.main, m.main)
# Load in our 2015 data --------------------------------------------------------
# First load in the GPS coordinates of our sites and subsites
gps_2015 <- read.csv(file = "raw_data/collection_sites_2015.csv")
# and load in the data on the subsites with the most butterflies,
# as we'll use those coordinates for the overall site. 
largest_subsite <- read.csv("processed_data/largest_subsite.csv")

# Keep only columns we want, for a left join later
# And get rid of the extra rows of subsites: use the 
# GPS coords of the subsite with the most butterflies
site_GPS <- gps_2015 %>%
  filter(subsite.name == "") %>%
  dplyr::select(site.collected = site.name, coord.N.decdeg, coord.W.decdeg)
subsite_GPS <- gps_2015 %>%
  filter(subsite.name %in% largest_subsite$subsite) %>%
  dplyr::select(site.collected = site.name, coord.N.decdeg, coord.W.decdeg)
subsite_GPS_melp <- gps_2015 %>%
  filter(subsite.name %in% largest_subsite_melp$subsite) %>%
  dplyr::select(site.collected = site.name, coord.N.decdeg, coord.W.decdeg)

# Combine and remove intermediates
gps_2015 <- rbind(site_GPS, subsite_GPS)
gps_2015_melp <- rbind(site_GPS, subsite_GPS_melp)
rm(largest_subsite, site_GPS, subsite_GPS)
rm(largest_subsite_melp, subsite_GPS_melp)

# First load in the individual data and sum up by site
# We use all the sites we collected at, 17 total. 
us <- read_csv(file = "raw_data/erato_collected_2015.csv") %>% 
  sumPhenoCounts(., "site.collected") %>%
  mutate(year = 2015) %>%
  left_join(., gps_2015, by = "site.collected") %>%
  dplyr::select(year, site.collected, coord.N.decdeg, coord.W.decdeg,
         A.melanized, B.hetero, C.west.col, D.postman)
rm(gps_2015)

# Combine, sort west-to-east, remove intermediates-------------------------

joint <- rbind(mallet, blum, us) %>%
  dplyr::select(year, site.collected, coord.N.decdeg, coord.W.decdeg,
                A.melanized, B.hetero, C.west.col, D.postman) %>%
  arrange(coord.W.decdeg)

rm(mallet, blum, us)


# Fit a cubic transect ----------------------------------------------------
# Fit the model, and save a summary of the model as a tidy df.
cubic.model <- lm(coord.N.decdeg ~ coord.W.decdeg + I(coord.W.decdeg^2) + I(coord.W.decdeg^3), data = joint)
cm.table <- tidy(lm(coord.N.decdeg ~ coord.W.decdeg + I(coord.W.decdeg^2) + I(coord.W.decdeg^3), data = joint))


# Find closest point on the transect --------------------------------------
# First, make a data table of x,y values that describe the fitted cline.
# This data frame will be used by the addPointonTransect function
cubic.longs <- data.frame(coord.W.decdeg = seq(-83, -77, 0.0002)) # eqivalent to 1 point every second
cubic.trans <- data.frame(cbind(coord.N.decdeg = predict(cubic.model, cubic.longs), cubic.longs))
rm(cubic.longs)

# Now, add the coordinates for the position along the cline
# This is a bit slow, but only needs to be done once, so little point in profiling.
joint_transect <- addPointOnTransect(joint) %>%
  arrange(tran.Coord.W)



# Add transect distance column --------------------------------------------
joint_transect <- addDistanceColsArc(joint_transect, cm.table$estimate)

write.csv(joint_transect, "processed_data/joint_transect.csv", row.names = F)

# Generate a set of points to use for forest analysis ---------------------
# I want to generate a latitude values along the transect such that
# each point is a set distance apart.
joint_transect <- read.csv("processed_data/joint_transect.csv")

# Then, I can generate a set of lat/long coordinates to use for the
# forest analysis. 

generated <- gen_points_on_trans(transect.coeffs = cubic.model$coefficients,
                                 start.lat =  joint_transect$tran.Coord.W[1],
                                 end.dist = 700, 
                                 interval = 15)

write.csv(generated, file = "processed_data/generated_sites_along_transect.csv", row.names = F)


