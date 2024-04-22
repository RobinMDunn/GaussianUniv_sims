# Simulate power of crossfit LRT.

suppressMessages(library(data.table))
suppressMessages(library(progress))
suppressMessages(library(MASS))

# Arguments for d (dimensional of Normal dist),
# number of observations, and number of simulations.
d_vec <- c(1, 2, 10)
n_obs <- 1000
n_sim <- 5000

# Set alpha level
alpha_val <- 0.1

# Construct data frame to store results
results <- data.table(n = n_obs, 
                      alpha = alpha_val,
                      expand.grid(d = d_vec, 
                                  theta_star = seq(0.001, 0.317, by = 0.001)),
                      power = NA_real_)

# Set up progress bar
pb <- progress_bar$new(format = "sim :current / :total [:bar] :eta",
                       total = nrow(results), clear = T, show_after = 0)

# Check whether each test point is inside usual LRT.
for(row in 1:nrow(results)) {
  
  # Increment progress bar
  pb$tick()
  
  # Extract n, d, and theta_star
  n_val <- results[row, n]
  
  d_val <- results[row, d]
  
  theta_star_val <- results[row, theta_star]
  
  theta_star_vec <- rep(theta_star_val, times = d_val)
  
  # Vector to indicate rejection of H_0
  reject_h0 <- rep(NA, n_sim)
  
  for(sim in 1:n_sim) {
    
    # Generate data
    Y <- matrix(mvrnorm(n = n_val, mu = theta_star_vec, Sigma = diag(d_val)),
                ncol = d_val)
    
    # Split Y into Y_0 and Y_1
    Y_0_indices <- sample(1:n_val, size = n_val / 2)
    
    Y_1_indices <- setdiff(1:n_val, Y_0_indices)
    
    Y_0 <- matrix(Y[Y_0_indices, ], ncol = d_val)
    
    Y_1 <- matrix(Y[Y_1_indices, ], ncol = d_val)
    
    # Overall means of Y_0 and Y_1 observations
    Y_0_bar <- apply(Y_0, MARGIN = 2, FUN = mean)
    
    Y_1_bar <- apply(Y_1, MARGIN = 2, FUN = mean)
    
    # Check whether you would reject H_0: theta_0 = 0
    reject_h0[sim] <- 
      as.numeric(.5*(exp(-n_val * sum((Y_0_bar - Y_1_bar)^2) / 4 +
                           n_val * sum((Y_0_bar)^2) / 4) +
                       exp(-n_val * sum((Y_0_bar - Y_1_bar)^2) / 4 +
                             n_val * sum((Y_1_bar)^2) / 4)) > 
                   1 / alpha_val)
  }
  
  # EStimated power of testing H_0: theta_0 = rep(0, d_val) at theta_star_vec
  results[row, power := mean(reject_h0)]
  
}

# Save simulation results.
fwrite(results, file = "sim_data/fig06_crossfit_LRT.csv")

