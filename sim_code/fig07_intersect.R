# Simulate power of confidence set intersection test.
# H_0: ||theta*|| \in [0.5, 1.0]
# H_1: ||theta*|| \notin [0.5, 1.0]

suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(progress))
suppressMessages(library(MASS))

# Read in arguments for d (dimensional of Normal dist),
# start coordinate, end coordinate, increment, 
# number of simulations, number of subsamples, and number of observations.

d_val <- 2
coord_start <- 0
coord_end <- 2
increment <- 0.1
n_sim <- 1000
n_obs <- 1000

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  args <- as.numeric(args)
  d_val <- args[1]
  coord_start <- args[2]
  coord_end <- args[3]
  increment <- args[4]
  n_sim <- args[5]
  n_obs <- args[6]
}

# Set alpha level
alpha_val <- 0.1

# Construct data frame to store results
results <- data.table(n = n_obs,
                      n_sim = n_sim,
                      alpha = alpha_val,
                      d = d_val,
                      expand.grid(sim = 1:n_sim,
                                  theta_star_1 = seq(coord_start, coord_end, 
                                                     by = increment)),
                      reject = NA_real_)

# Set up progress bar
pb <- progress_bar$new(format = "sim :current / :total [:bar] :eta",
                       total = nrow(results), clear = T, show_after = 0)

# Check whether each test point is inside usual LRT.
for(row in 1:nrow(results)) {
  
  # Increment progress bar
  pb$tick()
  
  # Extract sim and theta_star
  sim_val <- results[row, sim]
  
  theta_star_1_val <- results[row, theta_star_1]
  
  # Generate data
  Y <- cbind(matrix(rnorm(n = n_obs, mean = theta_star_1_val, sd = 1), 
                    ncol = 1, nrow = n_obs),
             matrix(rnorm(n = n_obs * (d_val - 1), mean = 0, sd = 1), 
                    ncol = d_val - 1, nrow = n_obs))
  
  # Get overall mean
  Y_bar <- apply(Y, MARGIN = 2, FUN = mean)
  
  # Get norm of Y_bar
  Y_bar_norm <- sqrt(sum(Y_bar^2))
  
  # Project mean onto S_1 \setdiff S_0.5
  if((Y_bar_norm >= 0.5) & (Y_bar_norm <= 1)) {
    theta <- Y_bar 
  } else if(Y_bar_norm > 1) {
    theta <- Y_bar / Y_bar_norm
  } else if(Y_bar_norm < 0.5) {
    theta <- 0.5 * Y_bar / Y_bar_norm
  }
  
  # theta is now the closest value to Y_bar inside S_1 \setdiff S_0.5.
  # Reject H_0 if theta \notin C_n^{LRT}(\alpha)
  
  cutoff <- qchisq(p = alpha_val, df = d_val, lower.tail = F) / n_obs
  
  results[row, reject := as.numeric(sum((theta - Y_bar)^2) > cutoff)]
  
}

# Save rejection proportions
results <- results %>% 
  group_by(n, n_sim, alpha, d, theta_star_1) %>% 
  dplyr::summarise(reject = mean(reject))

# Save simulation results.
fwrite(results, file = paste0("sim_data/fig07_intersect_", d_val,
                              "_", coord_start, ".csv"))

