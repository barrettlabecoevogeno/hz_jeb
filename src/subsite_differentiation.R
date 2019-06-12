############################################
## Tim Thurman                            ##
## Code for analysis of                   ##
## H. erato hybrid zone in eastern Panama ##
## Testing for differentiation between    ##
## subsites to ensure that they can be    ##
## merged together                        ##
############################################


# Load packages and functions----------------------------------------------
## Clear the environment and load in necessary packages
rm(list=ls())
library(tidyverse)
source("src/functions.R")


# Load in Data ------------------------------------------------------------
## Load in the raw, individual collection data.
individuals <- read_csv(file = "raw_data/erato_collected_2015.csv") 

# Get allele counts for each subsite-------------------------------------------------
allele_counts_subsite <- individuals %>%
  sumPhenoCounts(., c("subsite", "site.collected")) %>% # sum up # of each phenotype per subsite, keep the site column
  filter(is.na(subsite)==F) %>% # get rid of stuff that doesn't have a subsite
  phenoToAlleleCount() 

 
# Test for subsite differentiation ----------------------------------------
subsite_diff <- allele_counts_subsite %>%
  testSubsiteDiff()

# Results: In all sites where we did subsampling, no evidence that subsamples have different allele frequencies
# Can be merged.

# Tally subsite totals ----------------------------------------------------
# Later, will want to use the GPS coordinates of the subsite that has the most individuals
# for each site
largest_subsite <- data.frame(cbind(site.collected = unique(allele_counts_subsite$site.collected),
                                    subsite = rep(NA, times=length(unique(allele_counts_subsite$site.collected)))),
                              stringsAsFactors = F)

for (site in unique(allele_counts_subsite$site.collected)) {
  sub <- subset(allele_counts_subsite, site.collected == site)
  res.row <- which(largest_subsite$site.collected == site)
  max <- which(sub$totalAlleles == max(sub$totalAlleles))
  if (length(max) > 1) { # if there's a tie, randomly pick one subsite
    max <- sample(max, 1)
  }
  largest_subsite$subsite[res.row] <- sub$subsite[max]
}
write.csv(largest_subsite, file = "processed_data/largest_subsite.csv", row.names = F)
