
# Load pacakges -----------------------------------------------------------
rm(list=ls())
library(tidyverse)
library(broom)


# Load data ---------------------------------------------------------------
load("results/corr_simulation_res.Rdata")
load("results/bayes_simulation_res.Rdata")


# Analysis ----------------------------------------------------------------
full <- rbind(bayes_results, corr) %>% 
  arrange(model, data_set)
rm(bayes_results, corr)

simtime <- full %>% 
  group_by(model) %>% 
  summarise(mean.time = mean(time/60))

RMSD_by_parameter_set <- full %>% 
  mutate(difference = estimate - sim_val,
         sq.diff = (estimate - sim_val)^2) %>% 
  group_by(parameter, param_set, model) %>% 
  summarise(RMSD = sqrt(mean(sq.diff))) %>% 
  ungroup()

# And then do a t.test
pair.t.test <- RMSD_by_parameter_set %>% 
  spread(model, RMSD) %>% 
  group_by(parameter) %>% 
  do(tidy(t.test(.$bayes, .$corr.HZAR, paired = T))) %>% 
  ungroup() %>% 
  rename(df = parameter) %>% 
  mutate(parameter = c("center", "pmax", "pmin", "width"), 
         param_set = 0,
         compare = "average")

perr_correct_by_model <- full %>%
  mutate(greater.low = sim_val > low_est,
         less.high = sim_val < up_est) %>%
  mutate(within = (greater.low + less.high == 2)) %>%
  group_by(model) %>%
  summarize(correct = sum(within),
            total = length(within)) %>%
  ungroup() %>%
  mutate(perc.corr = correct/total)







