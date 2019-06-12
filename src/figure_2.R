############################################
## Tim Thurman                            ##
## Code for analysis of                   ##
## H. erato hybrid zone in eastern Panama ##
## Figure 2                               ##
############################################

# Load packages -----------------------------------------------------------

library(tidyverse)
library(bayesplot)
source("src/functions.R")
source("src/hzar_functions.R")

# Load and prep data ------------------------------------------------------
joint <- read.csv("processed_data/joint_transect.csv") %>% 
  mutate(AA = A.melanized, 
         Aa = B.hetero, 
         aa = C.west.col + D.postman,
         N = A.melanized + B.hetero + C.west.col + D.postman)

plot.dat <- joint %>%
  filter(year == "2015") %>% 
  mutate(alpha = 2*AA + Aa + 0.5,
         beta = (2*(AA + Aa + aa) - (2*AA + Aa)) + 0.5,
         obsP = (2*AA + Aa)/(2*(AA + Aa + aa))) %>%
  mutate(low95 = qbeta(0.025, alpha, beta),
         up95 = qbeta(0.975, alpha, beta))


# Using the cline equations to generate the best-fit lines for each year
# Blum was right tail
xrange <- 0:700  
blum.Best <- right.eqn.vec(xrange, 467.15, 59.52, 0.04, 0.95, 20.56, 0.57)

# Mallet was no tail
mallet.Best <- none.eqn.vec(xrange, 516.05, 52.87, 0.07, 0.91)


# Thurman was right tail
thurman.Best <- right.eqn.vec(xrange, 450.89, 93.18, 0.04, 0.94, 24.81, 0.66)

# Thurman
load("results/thurman_bestCline.Rdata")
thurman <- right
rm(right)


desired_params <- grep("\\[|_", names(thurman), invert = T, value = T)
post <- as.matrix(thurman, pars = desired_params)


# Make plot ---------------------------------------------------------------
a <- ggplot() +
  geom_blank(aes(x = transect.dist, y = obsP), data = plot.dat) +
  scale_x_continuous(breaks= seq(300, 600, 50)) +
  theme(axis.text = element_text(size = 10)) +
  xlab("Distance along transect (km)") +
  ylab(expression("Frequency of Cr"[HYD])) +
  theme_default() +
  theme(axis.text = element_text(size = 10))
for (samp in 1:post[1]) {
  a <- a + geom_line(aes_string(x = xrange, y = right.eqn.vec(xrange, post[samp, 1], post[samp, 2],
                                                              post[samp,3], post[samp, 4],
                                                              post[samp, 5], post[samp, 6])), color = alpha("black", 0.02))
}
b <- a + geom_linerange(size = 0.5, color = "black", aes(x = transect.dist, ymin = low95, ymax = up95), data = plot.dat) +
  geom_point(size = 2.5, color = "black", aes(x = transect.dist, y = obsP), data = plot.dat) +
  geom_line(aes_string(x = xrange, y = mallet.Best), color = "black", size = 1, linetype = "dotted") +
  geom_line(aes_string(x = xrange, y = blum.Best), color = "black", size = 1, linetype = "dashed") +
  geom_line(aes_string(x = xrange, y = thurman.Best), color = "#FF8C00", size = 1) +
  coord_cartesian(x = c(290, 620), y = c(0,1))
c <- b + 
  geom_text(aes_string(x = 600, y = 0.25, label = "1982"), size = 4) +
  geom_text(aes_string(x = 600, y = 0.18, label = "1999"), size = 4) +
  geom_text(aes_string(x = 600, y = 0.11, label = "2015"), size = 4) +
  geom_line(aes_string(x = c(556, 585), y = c(0.245, 0.245)), linetype = "dotted", size = 1) +
  geom_line(aes_string(x = c(556, 585), y = c(0.175, 0.175)), linetype = "dashed", size = 1) +
  geom_line(aes_string(x = c(556, 581), y = c(0.105, 0.105)), linetype = "solid", size = 1, color = "#FF8C00")
#c

ggsave("results/figure_2.pdf", c,  height = 90, width = 166, units = "mm", dpi = 1000)
