# Comparing results from HZAR
# to the new Bayesian method

# Load packages -----------------------------------------------------------
rm(list=ls())
library(tidyverse)
library(stringr)
library(hzar)
library(doMC)
registerDoMC()
library(rstan)
library(tictoc)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
source("src/functions.R")
source("src/hzar_functions.R")
correct_fis <- correct_Fis

# Simulate data in a list, so that it can be reused -----------------------

pos_fis <- c(0, 0.1, 0.25, 0.5, 0.75, 1)
nrep <- 15
center <- 200
width <- c(20, 80)
pmin <- c(0.04, 0.15)
pmax <- c(0.85, 0.97)

set.seed(92)
sim_dat_list <- list()
index <- 1
for (w in width) {
  for (pi in pmin) {
    for (px in pmax) {
      for (fis in pos_fis) {
        rep <- 1
        while (rep <= nrep) { # for the number of iterations
          simulated <- sim_data_from_cline(transect_distances = seq(0,400,10), n_ind = 40, Fis = fis,
                                           decrease = F, center = center, width = w, pmin = pi, pmax = px)
          sim_dat_list[[index]] <- list()
          sim_dat_list[[index]][[1]] <- as.character(center)
          sim_dat_list[[index]][[2]] <- as.character(w)
          sim_dat_list[[index]][[3]] <- as.character(pi)
          sim_dat_list[[index]][[4]] <- as.character(px)
          sim_dat_list[[index]][[5]] <- as.character(fis)
          sim_dat_list[[index]][[6]] <- as.character(rep)
          sim_dat_list[[index]][[7]] <- simulated
          index <- index + 1
          rep <- rep + 1
        }
      }
    }
  }
}

rm(simulated)

param_set_key <- expand.grid(width = width, pmin = pmin, pmax = pmax, sim_fis = pos_fis) %>% 
  arrange(width, pmin, pmax, sim_fis) %>% 
  mutate(param_set = 1:48)
data_set_key <- expand.grid(width = width, pmin = pmin, pmax = pmax, sim_fis = pos_fis, replicate = 1:15) %>% 
  left_join(param_set_key) %>% 
  arrange(param_set, replicate) %>% 
  mutate(data_set = 1:720)

# Fit bayesian clines for each element of the list ------------------------

results <- NULL
for (element in sim_dat_list[bayes_rerun4]) {
  c <- as.numeric(element[[1]])
  w <- as.numeric(element[[2]])
  pi <- as.numeric(element[[3]])
  px <- as.numeric(element[[4]])
  fis <- as.numeric(element[[5]])
  rep <- as.numeric(element[[6]])
  simulated <- element[[7]]
  
  print(paste("Fis = ", fis, " and replicate equals ", rep))
  
  # Fit the Bayesian model
  tic()
  none.stanfit <- stan("src/stan_models/multinomial/multi_free_none.stan",
                       data = list(N = dim(simulated)[1], 
                                   genos = as.matrix(simulated[,5:7]),
                                   transectDist = simulated$transectDist),
                       iter = 10000, warmup = 3000, chains = 4,
                       control = list(adapt_delta = 0.99),
                       init = list(list(center = abs(rnorm(n = 1, c, 20)), 
                                        width = abs(rnorm(1, w, 15)),
                                        pmin = runif(1, 0, 0.2),
                                        pmax = runif(1, 0.8, 1)),
                                   list(center = abs(rnorm(n = 1, c, 20)), 
                                        width = abs(rnorm(1, w, 15)),
                                        pmin = runif(1, 0, 0.2),
                                        pmax = runif(1, 0.8, 1)),
                                   list(center = abs(rnorm(n = 1, c, 20)), 
                                        width = abs(rnorm(1, w, 15)),
                                        pmin = runif(1, 0, 0.2),
                                        pmax = runif(1, 0.8, 1)),
                                   list(center = abs(rnorm(n = 1, c, 20)), 
                                        width = abs(rnorm(1, w, 15)),
                                        pmin = runif(1, 0, 0.2),
                                        pmax = runif(1, 0.8, 1))))
  time <- toc()
  elapsed <- time$toc - time$tic
  notes <- ""
  # Check for problems, add notes as necessary
  if (none.stanfit@sim$chains != 4) {
    notes <- str_glue(notes, "IMPROPER NUM OF CHAINS")
  }
  
  cline_sum <- cline_summary(none.stanfit)
 
  # Summarize it down to the form we want
  rep_bayes_res <- cline_sum %>% 
    mutate(center = c,
           width = w,
           pmin = pi,
           pmax = px,
           sim_fis = as.character(fis),
           replicate = as.character(rep),
           model = "bayes",
           parameter = as.character(param),
           sim_val = c(c, w, pi, px),
           NOTES = notes,
           time = elapsed) %>% 
    dplyr::select(center, width, pmin, pmax, sim_fis, 
                  replicate, model, parameter, sim_val, 
                  estimate = mean, low_est = low_0.95_HPDI, up_est = up_0.95_HPDI, 
                  n_eff, Rhat, NOTES, time)
  # add to the final results
  results <- rbind(results, rep_bayes_res)
  # clear memory space so results don't get repeated in the case of failed iterations
  rm(none.stanfit, cline_sum, rep_bayes_res, time)
}

bayes_results <- results
# these results are saved as results/bayes_simulation_res.Rdata


# Fit corrected HZAR to each element of the list ------------------------

results <- NULL
for (element in sim_dat_list) {
  c <- as.numeric(element[[1]])
  w <- as.numeric(element[[2]])
  pi <- as.numeric(element[[3]])
  px <- as.numeric(element[[4]])
  fis <- as.numeric(element[[5]])
  rep <- as.numeric(element[[6]])
  simulated <- element[[7]]
  
  print(paste("Fis = ", fis, " and replicate equals ", rep))
  
  corr_sim <- simulated %>% 
    mutate(transect.dist = transectDist,
           Cr.HYD = emp.p,
           Cr.nSamples = (N*2)/(1 + emp.f))
  
  attempt <- 1
  while (attempt <= 15) {
    tic() #start timer
    hzar.corr.fit <- fitNoneHZAR(corr_sim, c = c, w = w) # try the fit
    time <- toc() # end timer
    elapsed <- time$toc - time$tic
    if (elapsed < 30) {
      # If the fitting failed, increment attempt and retry
      attempt <- attempt + 1
      next
    } else {
      # If the fitting worked, get the results and move on
      hzar.corr.fit <- modelSelect(hzar.corr.fit)
      notes <- str_glue("took ", attempt, " attempts to fit")
      break
    }
  }
  
  # compile results
  corr_hzar_res <- extractModelParams(hzar.corr.fit) %>% 
    mutate(center = c,
           width = w,
           pmin = pi,
           pmax = px,
           sim_fis = as.character(fis),
           replicate = as.character(rep),
           model = "corr.HZAR",
           parameter = as.character(parameter),
           sim_val = c(c, w, pi, px),
           n_eff = NA,
           Rhat = NA,
           NOTES = notes,
           time = elapsed)  %>% 
    dplyr::select(center, width, pmin, pmax, sim_fis, 
                  replicate, model, parameter, sim_val, 
                  estimate, low_est = lower, up_est = upper, 
                  n_eff, Rhat, NOTES, time)
  
  results <- rbind(results, corr_hzar_res)
  results$parameter[which(results$parameter == "pMin")] <- "pmin"
  results$parameter[which(results$parameter == "pMax")] <- "pmax"
  rm(hzar.corr.fit, corr_hzar_res, time)
}
corr_results <- results
# these results are saved as results/corr_simulation_res.Rdata
