# Evaluate power of usual LRT.

suppressMessages(library(data.table))
suppressMessages(library(progress))
suppressMessages(library(MASS))

# Arguments for d (dimensional of Normal dist)
# and number of observations
d_vec <- c(1, 2, 10)
n_val <- 1000

# Set alpha level
alpha_val <- 0.1

# Construct data frame to store results
results <- data.table(alpha = alpha_val,
                      n = n_val,
                      expand.grid(d = d_vec, 
                                  theta_star = seq(0.001, 0.317, by = 0.001)),
                      power = NA_real_,
                      power_approx_2 = NA_real_)

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
  
  # What is the power of testing H_0: theta_0 = rep(0, d_val) at theta_star_vec
  power_val <- 
    pchisq(q = qchisq(p = alpha_val, df = d_val, lower.tail = F), 
           df = d_val, ncp = n_val * sum(theta_star_vec^2), lower.tail = F)
  
  results[row, power := power_val]
  
  # Second approximation to power, using Normal CDF
  c_alpha_d <- qchisq(p = alpha_val, df = d_val, lower.tail = F)
  
  theta_star_norm <- sqrt(sum(theta_star_vec^2))
  
  power_approx <- 
    pnorm(q = (d_val + n_val*theta_star_norm^2 - c_alpha_d) /
            sqrt(2 * (d_val + 2*n_val*theta_star_norm^2)),
          lower.tail = T)

  results[row, power_approx_2 := power_approx]
}

# Save simulation results.
fwrite(results, file = "data/power/02_normal_usual_LRT.csv")

