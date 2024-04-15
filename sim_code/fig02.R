# Get data on simulated and approximate analytical expressions for 
# subsampling as B goes to infinity.

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(MASS))
suppressMessages(library(progress))

# Read in number of observations and dimension
n <- 1000 
d <- 10   

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  args <- as.numeric(args)
  n <- args[1]
  d <- args[2]
}

# Simulate one dataset
set.seed(n + d)

Y <- mvrnorm(n, mu = rep(0, d), Sigma = diag(1, nrow = d, ncol = d))

Y_bar <- apply(Y, MARGIN = 2, FUN = mean)

# Find lower/upper c_val cutoffs at which test stat approx = 15. 
# Compute simulated and approx test stat at values within this range.
# This will show similarity of approx and simulations in 
# alpha = 0.1 region (test stat <= 10) and slightly beyond.
approx_test_stat_func <- function(c_val, n, d, Y_bar) {
  return(exp((3*n / 10) * sum((Y_bar - c_val)^2)) * (2/5)^(d/2))
}

lower_c_val <- 
  uniroot(f = function(x) approx_test_stat_func(x, n, d, Y_bar) - 15, 
          lower = -0.2, upper = 0, extendInt = "downX",
          tol = .Machine$double.eps)$root

upper_c_val <- 
  uniroot(f = function(x) approx_test_stat_func(x, n, d, Y_bar) - 15, 
          lower = 0, upper = 0.2, extendInt = "upX",
          tol = .Machine$double.eps)$root

# Create data table to store results from simulated average test statistics
results_sim <- 
  data.table(c_val = seq(lower_c_val, upper_c_val, length.out = 20),
             n = n, d = d, method = "Simulated\nexpectation", value = NA_real_)

# Set up progress bar
pb <- progress_bar$new(total = nrow(results_sim),
                       format = "Row :current / :total [:bar] :eta")

# Get estimate of expectation from averaging over simulations and
# analytical approximation to expectation, over a range of theta values.
for(row in 1:nrow(results_sim)) {
  
  pb$tick()
  
  theta <- c(rep(results_sim[row, c_val], d))
  
  sim <- rep(NA, 100000)
  
  # Get estimate of expectation through simulation
  for(i in 1:length(sim)){
    Y_1_indices <- sample(1:n, size = n / 2)
    
    Y_0_indices <- setdiff(1:n, Y_1_indices)
    
    Y_0 <- as.matrix(Y[Y_0_indices, ])
    
    Y_1 <- as.matrix(Y[Y_1_indices, ])
    
    Y_bar_0 <- apply(Y_0, MARGIN = 2, FUN = mean)
    Y_bar_1 <- apply(Y_1, MARGIN = 2, FUN = mean)
    
    sim[i] <- exp(-(n/4)*sum((Y_bar_0 - Y_bar_1)^2) + 
                    (n/4)*sum((Y_bar_0 - theta)^2))
    
  }
  
  results_sim[row, value := mean(sim)]
  
}

# Create data table with analytical approximations to test statistics,
# at more granulated c (or theta) values
results_analytical <- 
  data.table(c_val = seq(lower_c_val, upper_c_val, length.out = 200), 
             n = n, d = d, method = "Analytical\nexpression", 
             value = NA_real_)

results_analytical[, value :=
                     exp((3*n / 10) * sum((Y_bar - c_val)^2)) * (2/5)^(d/2),
                   by = seq_len(nrow(results_analytical))]

results <- rbind(results_sim, results_analytical)

#####################
##### Save data #####
#####################

fwrite(results, file = paste0("data/radius/08_ddim_limit_results_d_", d,
                              "_n_", n, ".csv"))
