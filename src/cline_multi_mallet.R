rm(list = ls())
library(rethinking)
library(tidyverse)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)



# Load the data -----------------------------------------------------------
# take in the genotype data, for now treating the west columbian allele
# as a demophoon allele
mallet.geno <- read.csv("processed_data/joint_transect.csv") %>% 
  filter(year == "1982") %>%
  mutate(AA = A.melanized, Aa = B.hetero, aa = C.west.col + D.postman)



# Run the multinomial models ----------------------------------------------
# Get seeds from random.org
# Use random draws from the prior distribution as the start value for each chain

# No tails ----------------------------------------------------------------
set.seed(437)
none <- stan("src/stan_models/multinomial/multi_free_none.stan",
             data = list(N = dim(mallet.geno)[1], 
                         genos = as.matrix(mallet.geno[,13:15]),
                         transectDist = mallet.geno$transect.dist),
             iter = 10000, warmup = 3000, chains = 4,
             control = list(adapt_delta = 0.99),
             init = list(list(center = rnorm(n = 1, 350, 150), 
                              width = rnorm(1, 50, 100),
                              pmin = runif(1, 0, 0.2),
                              pmax = runif(1, 0.8, 1)),
                         list(center = rnorm(n = 1, 350, 150), 
                              width = rnorm(1, 50, 100),
                              pmin = runif(1, 0, 0.2),
                              pmax = runif(1, 0.8, 1)),
                         list(center = rnorm(n = 1, 350, 150), 
                              width = rnorm(1, 50, 100),
                              pmin = runif(1, 0, 0.2),
                              pmax = runif(1, 0.8, 1)),
                         list(center = rnorm(n = 1, 350, 150), 
                              width = rnorm(1, 50, 100),
                              pmin = runif(1, 0, 0.2),
                              pmax = runif(1, 0.8, 1))))



# Left tail ---------------------------------------------------------------
set.seed(682)
left <- stan("src/stan_models/multinomial/multi_free_left.stan",
             data = list(N = dim(mallet.geno)[1], 
                         genos = as.matrix(mallet.geno[,13:15]),
                         transectDist = mallet.geno$transect.dist),
             iter = 10000, warmup = 3000, chains = 4,
             control = list(adapt_delta = 0.99),
             init = list(list(center = rnorm(n = 1, 350, 150), 
                              width = rnorm(1, 50, 100),
                              pmin = runif(1, 0, 0.2),
                              pmax = runif(1, 0.8, 1),
                              deltaL = rexp(1, 0.05),
                              tauL = runif(1,0,1)),
                         list(center = rnorm(n = 1, 350, 150), 
                              width = rnorm(1, 50, 100),
                              pmin = runif(1, 0, 0.2),
                              pmax = runif(1, 0.8, 1),
                              deltaL = rexp(1, 0.05),
                              tauL = runif(1,0,1)),
                         list(center = rnorm(n = 1, 350, 150), 
                              width = rnorm(1, 50, 100),
                              pmin = runif(1, 0, 0.2),
                              pmax = runif(1, 0.8, 1),
                              deltaL = rexp(1, 0.05),
                              tauL = runif(1,0,1)),
                         list(center = rnorm(n = 1, 350, 150), 
                              width = rnorm(1, 50, 100),
                              pmin = runif(1, 0, 0.2),
                              pmax = runif(1, 0.8, 1),
                              deltaL = rexp(1, 0.05),
                              tauL = runif(1,0,1))))

# Right tail --------------------------------------------------------------
set.seed(825)
right <- stan("src/stan_models/multinomial/multi_free_right.stan",
              data = list(N = dim(mallet.geno)[1], 
                          genos = as.matrix(mallet.geno[,13:15]),
                          transectDist = mallet.geno$transect.dist),
              iter = 10000, warmup = 3000, chains = 4,
              control = list(adapt_delta = 0.99),
              init = list(list(center = rnorm(n = 1, 350, 150), 
                               width = rnorm(1, 50, 100),
                               pmin = runif(1, 0, 0.2),
                               pmax = runif(1, 0.8, 1),
                               deltaR = rexp(1, 0.05),
                               tauR = runif(1,0,1)),
                          list(center = rnorm(n = 1, 350, 150), 
                               width = rnorm(1, 50, 100),
                               pmin = runif(1, 0, 0.2),
                               pmax = runif(1, 0.8, 1),
                               deltaR = rexp(1, 0.05),
                               tauR = runif(1,0,1)),
                          list(center = rnorm(n = 1, 350, 150), 
                               width = rnorm(1, 50, 100),
                               pmin = runif(1, 0, 0.2),
                               pmax = runif(1, 0.8, 1),
                               deltaR = rexp(1, 0.05),
                               tauR = runif(1,0,1)),
                          list(center = rnorm(n = 1, 350, 150), 
                               width = rnorm(1, 50, 100),
                               pmin = runif(1, 0, 0.2),
                               pmax = runif(1, 0.8, 1),
                               deltaR = rexp(1, 0.05),
                               tauR = runif(1,0,1))))

# Mirror tail -------------------------------------------------------------
set.seed(322)
mirror <- stan("src/stan_models/multinomial/multi_free_mirror.stan",
               data = list(N = dim(mallet.geno)[1], 
                           genos = as.matrix(mallet.geno[,13:15]),
                           transectDist = mallet.geno$transect.dist),
               iter = 10000, warmup = 3000, chains = 4,
               control = list(adapt_delta = 0.99),
               init = list(list(center = rnorm(n = 1, 350, 150), 
                                width = rnorm(1, 50, 100),
                                pmin = runif(1, 0, 0.2),
                                pmax = runif(1, 0.8, 1),
                                deltaM = rexp(1, 0.05),
                                tauM = runif(1,0,1)),
                           list(center = rnorm(n = 1, 350, 150), 
                                width = rnorm(1, 50, 100),
                                pmin = runif(1, 0, 0.2),
                                pmax = runif(1, 0.8, 1),
                                deltaM = rexp(1, 0.05),
                                tauM = runif(1,0,1)),
                           list(center = rnorm(n = 1, 350, 150), 
                                width = rnorm(1, 50, 100),
                                pmin = runif(1, 0, 0.2),
                                pmax = runif(1, 0.8, 1),
                                deltaM = rexp(1, 0.05),
                                tauM = runif(1,0,1)),
                           list(center = rnorm(n = 1, 350, 150), 
                                width = rnorm(1, 50, 100),
                                pmin = runif(1, 0, 0.2),
                                pmax = runif(1, 0.8, 1),
                                deltaM = rexp(1, 0.05),
                                tauM = runif(1,0,1))))


# Ind tails ---------------------------------------------------------------
set.seed(854)
ind <- stan("src/stan_models/multinomial/multi_free_ind.stan",
            data = list(N = dim(mallet.geno)[1], 
                        genos = as.matrix(mallet.geno[,13:15]),
                        transectDist = mallet.geno$transect.dist),
            iter = 10000, warmup = 3000, chains = 4,
            control = list(adapt_delta = 0.99),
            init = list(list(center = rnorm(n = 1, 350, 150), 
                             width = rnorm(1, 50, 100),
                             pmin = runif(1, 0, 0.2),
                             pmax = runif(1, 0.8, 1),
                             deltaL = rexp(1, 0.05),
                             tauL = runif(1,0,1),
                             deltaR = rexp(1, 0.05),
                             tauR = runif(1,0,1)),
                        list(center = rnorm(n = 1, 350, 150), 
                             width = rnorm(1, 50, 100),
                             pmin = runif(1, 0, 0.2),
                             pmax = runif(1, 0.8, 1),
                             deltaL = rexp(1, 0.05),
                             tauL = runif(1,0,1),
                             deltaR = rexp(1, 0.05),
                             tauR = runif(1,0,1)),
                        list(center = rnorm(n = 1, 350, 150), 
                             width = rnorm(1, 50, 100),
                             pmin = runif(1, 0, 0.2),
                             pmax = runif(1, 0.8, 1),
                             deltaL = rexp(1, 0.05),
                             tauL = runif(1,0,1),
                             deltaR = rexp(1, 0.05),
                             tauR = runif(1,0,1)),
                        list(center = rnorm(n = 1, 350, 150), 
                             width = rnorm(1, 50, 100),
                             pmin = runif(1, 0, 0.2),
                             pmax = runif(1, 0.8, 1),
                             deltaL = rexp(1, 0.05),
                             tauL = runif(1,0,1),
                             deltaR = rexp(1, 0.05),
                             tauR = runif(1,0,1))))


# Look at results ---------------------------------------------------------
precis(none)
precis(left)
precis(right)
precis(mirror)
precis(ind)

# Compare parameters across models
coeftab(none, left, right, mirror, ind)

# Use WAIC to compare model fit within multinomial
compare(none, left, right, mirror, ind)


# Save the best fit model for future use ----------------------------------
# Here, the model with tht highest weight is the one with no introgression tails.
save(none, file = "results/mallet_bestCline.Rdata")

# and save all model results for other uses
save(none, left, right, mirror, ind, file = "results/mallet_allClines.Rdata")

