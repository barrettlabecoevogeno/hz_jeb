## Examining changes in cline position and width through time

# Load in packages and source functions -----------------------------------
rm(list=ls())
library(tidyverse)
library(rethinking)
library(bayesplot)
library(rstan)
library(stringr)
library(assertthat)
library(ggrepel)
library(coda)
source("src/functions.R")

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


# Get the results in number form, generate a table ------------------------

make_coef_tbl <- function(stanfit) {
  desired_params <- grep("\\[|_", names(stanfit), invert = T, value = T)
  post <- as.data.frame(stanfit) %>%
    dplyr::select(desired_params)
  params <- colnames(post)
  means <- post %>%
    summarise_all(mean) %>%
    t(.) %>%
    as.vector(.)
  ci <- as.matrix(post) %>%
    apply(X = ., MARGIN = 2, FUN = HPDI, prob = .95) %>%
    t(.)
  results <- data.frame(param = params, mean = round(means, digits = 2),
                        lower95 = round(ci[,1], digits = 2), upper95 = round(ci[,2], digits = 2))
  rownames(results) <- NULL
  return(results)
}

make_coef_tbl(mallet)
make_coef_tbl(blum)
make_coef_tbl(thurman)


# Cline movement ----------------------------------------------------------
center.m.post <- mallet %>% 
  as.data.frame(.) %>% 
  dplyr::select(center) %>% 
  as.matrix(.)
center.b.post <- blum %>% 
  as.data.frame(.) %>% 
  dplyr::select(center) %>% 
  as.matrix(.)
center.t.post <- thurman %>% 
  as.data.frame(.) %>% 
  dplyr::select(center) %>% 
  as.matrix(.)

post.shift.mallet.blum <- center.b.post - center.m.post
post.shift.blum.thurman <- center.t.post - center.b.post

mean(post.shift.mallet.blum)
mean(post.shift.blum.thurman)
HPDinterval(as.mcmc(post.shift.mallet.blum))
HPDinterval(as.mcmc(post.shift.blum.thurman))


# Change in width ---------------------------------------------------------
width.m.post <- mallet %>% 
  as.data.frame(.) %>% 
  dplyr::select(width) %>% 
  as.matrix(.)
width.b.post <- blum %>% 
  as.data.frame(.) %>% 
  dplyr::select(width) %>% 
  as.matrix(.)
width.t.post <- thurman %>% 
  as.data.frame(.) %>% 
  dplyr::select(width) %>% 
  as.matrix(.)

post.change.width.mallet.blum <- width.b.post - width.m.post
post.change.width.blum.thurman <- width.t.post - width.b.post
post.change.width.mallet.thurman <- width.t.post - width.m.post

mean(post.change.width.mallet.blum)
mean(post.change.width.blum.thurman)
mean(post.change.width.mallet.thurman)
HPDinterval(as.mcmc(post.change.width.mallet.blum))
HPDinterval(as.mcmc(post.change.width.blum.thurman))
HPDinterval(as.mcmc(post.change.width.mallet.thurman))
sum(post.change.width.blum.thurman > 0)/length(post.change.width.blum.thurman)



# Movement, width and selection -------------------------------------------

# CLEAN UP THIS SECTION. 
# From Blum 2002, the width of the cline should be:

# We can do some math from equations (2) and (3) of Blum 2002
# and derive equations which relate dispersal distance and selection
# to cline width and movement.

# Specifically, if s is selection, v is velocity per generation (assume 4 per year)
# w is cline width, and disp is dispersal distance:
calc_s <- function(v, w) {
  (20*v)/w
}
# if 16 in eqn 4 50
calc_s <- function(v, w) {
  (50*v)/w
}



# if 8 in eqn 4, 0.16
calc_disp <- function(v, w) {
  (((v^2)*(w^2))/.16)^.25
}
# if 16 in eqn 4, 0.32
calc_disp <- function(v, w) {
  (((v^2)*(w^2))/.32)^.25
}






# if 16 in eqn 4, 0.32
calc_s <- function(v, w) {
  sqrt(9600)*(v/w)
}
calc_disp <- function(v, w) {
  (((v^2)*(w^2))/3.84)^.25
}



# We can use our point estimates of movement 
# and width to estimate s and disp
v <- 0.25
w <- 93

s <- calc_s(v,w)
est_disp <- calc_disp(v,w)


width <- 4*sqrt((12*est_disp^2)/s)
width == w
velocity <- 0.1*sqrt(2*disp^2*s)
velocity == v


# Could thus use the full posterior of our width estimate, and full posterior of our 
# velocity estimate, to get the distribution of these parameters. 
# the velocity per year is the overall movement (shift from blum to thurman) divided by # of years
# and then 4 per y4
v_post <- abs(post.shift.blum.thurman/16/4)
w_post <- as.matrix(thurman, pars = "width")/4
est_s <- calc_s(v_post, w_post)
est_disp <- calc_disp(v_post, w_post)

width <- 4*sqrt((12*est_disp^2)/est_s)
velocity <- 0.1*sqrt(2*est_disp^2*est_s)*16*4

unique(round(width, digits = 2) == round(w_post, digits = 2))
unique(round(velocity, digits = 4) == round(velocity, digits = 4))

mean(est_s)
mean(est_disp)
HPDinterval(as.mcmc(est_s))
HPDinterval(as.mcmc(est_disp))

v_post_early <- abs(post.shift.mallet.blum/17/4)
w_post_early <- as.matrix(blum, pars = "width")/4
est_s_early <- calc_s(v_post_early, w_post_early)
est_disp_early <- calc_disp(v_post_early, w_post_early)

mean(est_s_early)
mean(est_disp_early)
HPDinterval(as.mcmc(est_s_early))
HPDinterval(as.mcmc(est_disp_early))

library(bahz)

# after a mini-scare, the equation from blum is correct,
# and my width is still 1/max slope. 
joint %>% 
  filter(year == 2015) %>% 
  rename(transectDist = transect.dist) %>% 
  plot_geno_cline(stanfit = thurman, data = .)
abline(a = -4.329, b = 1/93.18)
abline(a = -18.817, b = 4/93.18)

