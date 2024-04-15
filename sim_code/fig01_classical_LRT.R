# Simulations based on 
# "Comparing LRT to split LRT for the multivariate Normal case" section of 
# "Universal Inference Using the Split Likelihood Ratio Test."
# Evaluate coverage of usual LRT over a grid to compare performance.

suppressMessages(library(data.table))
suppressMessages(library(progress))
suppressMessages(library(MASS))

# Read in arguments for d (dimensional of Normal dist),
# start/end sim, and number of observations
d_val <- 2
start_sim <- 1
end_sim <- 6
n_obs <- 1000

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  args <- as.numeric(args)
  start_sim <- args[1]
  end_sim <- args[2]
  n_obs <- args[3]
}

# Set alpha level
alpha <- 0.1

# Set theta vectors
theta <- c(0, 0)

# Construct data frame to store results
results <- data.table(n = n_obs, 
                      d = d_val,
                      coord1 = NA_real_,
                      coord2 = NA_real_,
                      sim = rep(1:6, each = 501 * 501),
                      theta_star_coord1 = rep(theta[1], each = 501 * 501),
                      theta_star_coord2 = rep(theta[2], each = 501 * 501),
                      covered = NA_real_)

# Extract subset of results for given simulation
results <- results[start_sim <= sim & sim <= end_sim]

# Create grid based on given theta values
for(i in start_sim:end_sim) {
  
  theta_star_coord1_val <- unique(results[sim == i, theta_star_coord1])
  
  theta_star_coord2_val <- unique(results[sim == i, theta_star_coord2])
  
  results[sim == i, coord1 := rep(seq(theta_star_coord1_val - 0.25, 
                                  theta_star_coord1_val + 0.25, by = 0.001), 
                                  each = 501)]
  
  results[sim == i, coord2 := rep(seq(theta_star_coord2_val - 0.25, 
                                      theta_star_coord2_val + 0.25, by = 0.001), 
                                  times = 501)]
}

# Set up progress bar
pb <- progress_bar$new(format = "sim :current / :total [:bar] :eta",
                       total = nrow(results[sim == 1]), clear = T, show_after = 0)

# Generate n d-dimensional vectors Y. (Each observation is a row.)
set.seed(20200324)
Y <- mvrnorm(n = n_obs, mu = theta, Sigma = diag(d_val))

# Overall mean of Y observations
Y_bar <- apply(Y, MARGIN = 2, FUN = mean)

# Check whether each test point is inside usual LRT.
for(row in 1:nrow(results[sim == 1])) {
  
  # Increment progress bar
  pb$tick()
  
  # Extract n, d, sim, and theta
  n_val <- results[row, n]
  
  d_val <- results[row, d]
  
  sim_val <- results[row, sim]
  
  # Extract coordinates for test point on grid
  test_point <- as.numeric(results[row, .(coord1, coord2)])
  
  # Does usual LRT universal confidence set cover true theta?
  results[row, covered := sum((test_point - Y_bar)^2) <= 
            qchisq(p = 1 - alpha, df = d_val) / n_val]

}

# Repeat sim 1 results for sims 2-6.
results[sim == 2, covered := results[sim == 1, covered]]
results[sim == 3, covered := results[sim == 1, covered]]
results[sim == 4, covered := results[sim == 1, covered]]
results[sim == 5, covered := results[sim == 1, covered]]
results[sim == 6, covered := results[sim == 1, covered]]

# Save simulation results.
fwrite(results, 
       file = paste0("data/region/fixed_theta_fixed_data/",
                     "00_normal_usual_LRT_sim_",
                     start_sim, "_", end_sim, ".csv"))

