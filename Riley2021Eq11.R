# Equation 11 in Riley et al. (2021) (). Minimum sample size for external validation of a clinical prediction model with a binary outcome. Statistics in Medicine 40:4230–4251. (https://onlinelibrary.wiley.com/doi/full/10.1002/sim.9025) 
# We provide a solution to re-arranging the equation to give N as a function of SE(C), C and phi.

###---- Libraries ----###
library(dplyr)      # data wrangling
library(tidyr)      # for creating a grid of input combinations
library(magrittr)   # piping operator(s)
library(purrr)      # for list mapping functions
library(ggplot2)    # plotting
library(viridis)    # color schemes
library(patchwork)  # plot combination and layout 



eq.N <- function(C, phi, SE) {
  # Estimate sample size as a function of the C-statistic, SE(C) and phi (the outcome proportion)
  # Args:
  #  C = the C-statistic
  #  phi = the outcome proportion
  #  SE = the the standard error (SE) of the C-statistic
  # Returns:
  #  an estimate of the sample size, conditional on C, SE and phi
  #
  #  alpha, beta, mu and SE2 are calculated internally by the function. These are just intermediate steps:
  SE2 <- SE^2
  alpha <- C - 1
  beta <- 16*(-2 + C - 2*C^3 + C^4)*SE2*phi
  mu <- C*(-1 + C)
  
  numerator <- (C + C^2 - 4*C^3 + 2*C^4 + sqrt(mu*(beta + (1 - 2*mu)^2*mu - beta*phi)))
  denominator <- (4*SE2*(-1 + alpha)*(2 + alpha)*(-1 + phi)*phi)
  
  n <- numerator/denominator
  
  # If desired, n can be conservatively rounded up to the next integer:
  # n <- ceiling(n)
  
  return(n)
}

# Examples of use: 
# NOTE: in the paper, they report use of SE = 0.0255, but the Stata simulation code in the Supplementary material uses 0.02551
# 1st example from Section 3.3.1 of the paper: they report 1154
eq.N(C = 0.7, phi = 0.1, SE = 0.02551) %>% ceiling()  # 1154
# 2nd example from Section 3.3.1 of the paper: they report 302
eq.N(C = 0.8, phi = 0.5, SE = 0.02551) %>% ceiling()  # 302



###---- Section 4.3 of the paper ----###
# C = 0.8, SE = 0.0255, phi = 0.018; they reported estimated n = 4252
eq.N(C = 0.8, phi = 0.018, SE = 0.02551) %>% ceiling()     # n = 4252

# For the sensitivity analysis they used:
# (i) C = 0.75, SE = 0.0255, phi = 0.018; they reported estimated n = 5125
eq.N(C = 0.75, phi = 0.018, SE = 0.02551) %>% ceiling()         # n = 5125

# (ii) C = 0.85, SE = 0.0255, phi = 0.018; they reported estimated n = 3271
eq.N(C = 0.85, phi = 0.018, SE = 0.02551) %>% ceiling()          # n = 3271



###---- Reproduce the plots (panels) shown in Fig. 2 of the paper ----###
# Prep the input data:
SE <- seq(0.005, 0.085, 0.001) 
C <- seq(0.6, 0.9, 0.1)
phi <- seq(0.1, 0.5, 0.1)

# The input data frame:
x <- 
  # All combinations of SE, phi and C:
  tidyr::crossing(SE, phi, C) %>% 
  # sample size estimation:
  dplyr::mutate(n = purrr::pmap_dbl(list(C, phi, SE), eq.N)) %>% 
  # Set phi as a factor for plotting:
  dplyr::mutate(phi = factor(phi))


# As the panels of Fig. 2 vary only in the input value for C, we wrap the plotting code in a function:
plot.fxn <- function(myC) {
  # Plot SE(C) as a function of sample size for given inputs of the C-statistic
  # Args:
  #  myC = numeric value of the C-statistic
  # Returns:
  #  a ggplot graphic
  #
  x %>% 
    # filter so that the x-axis range is about the same as shown in Fig. 2:
    # notice the curly braces around myC so that it is passed into the filtering construct: 
    dplyr::filter(C == {{myC}}, n <= 1500) %>% 
    ggplot(., aes(x = n, y = SE, color = phi, group = phi)) +
    geom_line() +
    # use the same x-axis labels as in Fig. 2:
    scale_x_continuous(breaks = seq(200, 1400, 200)) +
    scale_color_viridis(discrete = TRUE) +
    theme_light() + 
    labs(x = "sample size", y = "SE(C)") +
    geom_hline(yintercept = 0.0255, linetype = "solid", color = "grey50") +
    # add the plot title:
    ggtitle(paste("c-statistic =", myC)) +
    theme(
      # remove the minor grid lines for the x axis:
      panel.grid.minor.x = element_blank(),
      # center the plot title:
      plot.title = element_text(hjust = 0.5),
      # place the legend inside the plot:
      legend.position = "inside",
      legend.position.inside = c(.95, .95),  
      legend.justification = c("right", "top"))
}

# Generate the plots
p1 <- plot.fxn(myC = 0.6)
p2 <- plot.fxn(myC = 0.7)
p3 <- plot.fxn(myC = 0.8)
p4 <- plot.fxn(myC = 0.9)


# patchwork is used to lay out the plots in a 2 by 2 grid, place a common legend for phi at the bottom
(p1 + p2)/(p3 + p4) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
