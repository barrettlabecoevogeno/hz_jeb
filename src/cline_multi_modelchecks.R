## Doing model diagnostics and
## Posterior predictive checks
## For the best-fit models from each year

# Load in packages and source functions -----------------------------------
rm(list=ls())
library(tidyverse)
library(bayesplot)
library(stringr)
library(coda)
library(lemon)

# Load in best-fit models -------------------------------------------------
# Mallet
load("results/mallet_bestCline.Rdata")
mallet <- none
rm(none)

# Blum
load("results/blum_bestCline.Rdata")
blum <- right
rm(right)

# Thurman
load("results/thurman_bestCline.Rdata")
thurman <- right
rm(right)


# Load in Joint transect --------------------------------------------------
joint <- read.csv("processed_data/joint_transect.csv") %>% 
  mutate(AA = A.melanized, Aa = B.hetero, aa = C.west.col + D.postman)

# Model checking: convergence with Rhat -----------------------------------
# Checking to see if the Rhat convergence diagnostic deviates from 1 for any parameters
rhat_all <- c(rhat(mallet), rhat(blum), rhat(thurman))
which(rhat_all > 1.01) # All fine, nothing above 1.01

# Model checking: effective samples ---------------------------------------
# Looking at the ratio of efective samples to total samples.
neff_ratio_all <- c(neff_ratio(mallet), neff_ratio(blum), neff_ratio(thurman))
which(neff_ratio_all < 0.1) # nothing less than 0.1, good. 

# Model checking: trace plots ---------------------------------------------

# Visually examine the trace plots for the parameters we care about
# Don't plot the traces of expected p, y_rep, deviance, and logliklihood


plot_trace <- function(stanfit) { # Take in a stanfit model
  # convert it to an array, which Bayesplot likes for plotting
  post <- as.array(stanfit)
  # extract the NUTS parameters for diagnosos
  np <- nuts_params(stanfit) 
  # Get the name of the stanfit model this is being done on
  name <- str_to_title(str_replace(deparse(substitute(stanfit)), ".eff", ""))
  # Use a high-contrast color scheme
  color_scheme_set("viridis")
  # Then plot the traces using mcmc_trace, all in one column
  keep <- grep("\\[|_", names(stanfit), invert = T, value = T)
  mcmc_trace(post, pars = keep, np = np,
             facet_args = list(ncol = 1, strip.position = "left")) +
    ggtitle(label = paste(name, "cline, model=", stanfit@model_name, sep = " "))
}

plot_trace(mallet)
plot_trace(blum)
plot_trace(thurman)

# Posterior predictive checks ---------------------------------------------
postPred_intervals <- function(stanfit, raw.data.frame, prob) {
  ppc_intervals <- as.data.frame(stanfit) %>%
    dplyr::select(starts_with("y_rep")) %>%
    as.mcmc(.) %>%
    HPDinterval(., prob = prob) %>% 
    as.data.frame(.) %>% 
    mutate(param = row.names(.)) %>% 
    separate(param, into = c("y", "rep", "site", "genotype")) %>% 
    dplyr::select(site, genotype, ppc_low = lower, ppc_up = upper)
  ppc_intervals$genotype[which(ppc_intervals$genotype == 1)] <- "AA"
  ppc_intervals$genotype[which(ppc_intervals$genotype == 2)] <- "Aa"
  ppc_intervals$genotype[which(ppc_intervals$genotype == 3)] <- "aa"

  raw.data.frame %>% 
    dplyr::select(transect.dist, AA, Aa, aa, site.collected) %>% 
    dplyr::mutate(total = AA + Aa + aa,
                  site = as.character(seq(from = 1, to = dim(.)[1], by = 1))) %>% 
    gather(AA:aa, key = "genotype", value = "obs.count") %>% 
    left_join(ppc_intervals, by = c("site", "genotype")) %>% 
    mutate(inRange = ifelse(obs.count <= ppc_up & obs.count >= ppc_low, T, F))
}

mallet.ppc <- postPred_intervals(mallet, filter(joint, year == "1982"), .95)
blum.ppc <- postPred_intervals(blum, filter(joint, year == "1999"), .95)
thurman.ppc <- postPred_intervals(thurman, filter(joint, year == "2015"), .95)



ppc_res <- mallet.ppc
plot_ppc <- function(ppc_res, ...) {
  cols <- c("FALSE" = "darkorange", "TRUE" = rgb(122,207,221, maxColorValue = 255))
  ppc_res <- droplevels(ppc_res)
  
  ppc_res$site.collected <- factor(ppc_res$site.collected, levels = unique(ppc_res$site.collected))
  ppc_res$genotype <- factor(ppc_res$genotype, levels = c("AA", "Aa", "aa"))
  ppc_res %>% 
    mutate(obs = obs.count/total,
           ppc_low = ppc_low/total,
           ppc_up = ppc_up/total) %>% 
    ggplot(aes(x = site.collected, ymin = ppc_low, ymax = ppc_up, y = obs, color = inRange)) +
    geom_linerange(position = position_dodge(width = 1)) +
    geom_point() + 
    scale_color_manual(values = cols) +
    facet_rep_grid(genotype~.) +
    theme_default() +
    theme(legend.position="none",
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.45)) + 
    xlab("Site") +
    ylab("Genotype frequency")
}

plot_ppc(mallet.ppc)
plot_ppc(blum.ppc)
plot_ppc(thurman.ppc)
