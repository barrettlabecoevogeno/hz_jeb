############################################
## Tim Thurman                            ##
## Code for analysis of                   ##
## H. erato hybrid zone in eastern Panama ##
## Calculating allele frequencies for     ## 
## the three alleles                      ##
############################################



# Load packages and functions----------------------------------------------
## Clear the environment and load in necessary packages
rm(list=ls())
library(tidyverse)
library(readxl)
library(assertthat)
source("src/functions.R")


# Load in joint transect --------------------------------------------------

joint <- read.csv("processed_data/joint_transect.csv") 


# Frequencies of all 3 alleles --------------------------------------------
# From Mallet 1985 appendix, the ratio of the frequencies of the Central american allele to the 
# west columbian allele (ie., p(CA)/p(WC)) is equal to
# (f(D) + f(D)*(f(d) + f(c)))/f(c), 
# where f(D) is the frequncy of type D (full postman) individuals and f(C) is the frequency of West Colombian individuals
# These are maximum likelihood estimates of each allele frequency
# Only relevant for sites where there was a west colombian homozygote
# See  yel_part_proof in the docs folder for a more detailed explanation of the method,
# Which is just a simple algebraic rearrangement. 

# When both yellow homozygotes are present, do the full thing to partition
both_yellow_present <- joint %>%
  filter(C.west.col > 0) %>%
  filter(D.postman > 0) %>%
  mutate(n.Inds = A.melanized + B.hetero + C.west.col + D.postman) %>%
  mutate(freq.C = C.west.col/n.Inds,
         freq.D = D.postman/n.Inds) %>%
  mutate(freqHYDallele = (2*A.melanized + B.hetero)/(2*n.Inds)) %>%
  mutate(z = freq.D + freq.D*(freq.C + freq.D)) %>%
  mutate(freqCAallele = (z*(1-freqHYDallele))/(z + freq.C)) %>%
  mutate(freqWCallele = 1 - freqHYDallele - freqCAallele) %>%
  dplyr::select(site.collected, freqHYDallele, freqCAallele, freqWCallele, year, transect.dist)

# When west colombian yellow types are present, but no CA,
# assume all yellow alleles are from west colombia
only_westcol <- joint %>% 
  filter(C.west.col > 0) %>%
  filter(D.postman == 0) %>%
  mutate(n.Inds = A.melanized + B.hetero + C.west.col + D.postman)%>%
  mutate(freqHYDallele = (2*A.melanized + B.hetero)/(2*n.Inds)) %>%
  mutate(freqCAallele = 0) %>%
  mutate(freqWCallele = 1 - freqHYDallele) %>%
  dplyr::select(site.collected, freqHYDallele, freqCAallele, freqWCallele, year, transect.dist)

# And the opposite: when CA yellow types are present, but no WC,
# assume all yellow alleles are from  Central america
only_ca <- joint %>%
  filter(C.west.col == 0) %>%
  filter(D.postman > 0) %>%
  mutate(n.Inds = A.melanized + B.hetero + C.west.col + D.postman) %>%
  mutate(freqHYDallele = (2*A.melanized + B.hetero)/(2*n.Inds)) %>%
  mutate(freqCAallele = 1 - freqHYDallele) %>%
  mutate(freqWCallele = 0) %>%
  dplyr::select(site.collected, freqHYDallele, freqCAallele, freqWCallele, year, transect.dist)
  
# When neither yellow homozygote is present, assume any yellow alleles (in hets)
# are from central America: it is the more common allele overall
no_yel_homo <- joint %>%
  filter(C.west.col == 0) %>%
  filter(D.postman == 0) %>%
  mutate(n.Inds = A.melanized + B.hetero + C.west.col + D.postman) %>%
  mutate(freqHYDallele = (2*A.melanized + B.hetero)/(2*n.Inds)) %>%
  mutate(freqCAallele = 1 - freqHYDallele) %>%
  mutate(freqWCallele = 0) %>%
  dplyr::select(site.collected, freqHYDallele, freqCAallele, freqWCallele, year, transect.dist)


## Combine the 4 situations, calculate allele frequencies, round to 2 digits
triAllelic_freqs <- rbind(both_yellow_present, no_yel_homo, only_ca, only_westcol) %>%
  mutate(f.HYD = round(freqHYDallele, digits = 2)) %>%
  mutate(f.CA = round(freqCAallele, digits = 2)) %>%
  mutate(f.WC = round(freqWCallele, digits = 2)) %>%
  dplyr::select(site.collected, f.HYD, f.CA, f.WC, year, transect.dist)

# remove intermediates
rm(both_yellow_present, no_yel_homo, only_ca, only_westcol)

write.csv(triAllelic_freqs, "results/triallelic_freqs.csv", row.names = F)


