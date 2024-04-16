# Simulations based on 
# "Comparing LRT to split LRT for the multivariate Normal case" section of 
# "Universal Inference Using the Split Likelihood Ratio Test."
# Evaluate coverage of crossfit LRT over a grid.

suppressMessages(library(data.table))
suppressMessages(library(progress))
suppressMessages(library(MASS))

# Read in arguments for d (dimensional of Normal dist),
# sim number, and number of observations
d_val <- 2
sim_val <- 1
n_val <- 1000

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  args <- as.numeric(args)
  sim_val <- args[1]
  n_val <- args[2]
}

# Set alpha level
alpha <- 0.1

# Set theta vectors
theta <- c(0, 0)

# Construct data frame to store results
results <- data.table(n = n_val, 
                      d = d_val,
                      coord1 = NA_real_,
                      coord2 = NA_real_,
                      sim = rep(sim_val, each = 501 * 501),
                      Ybar0_coord1 = NA_real_,
                      Ybar0_coord2 = NA_real_,
                      Ybar1_coord1 = NA_real_,
                      Ybar1_coord2 = NA_real_,
                      theta_star_coord1 = theta[1],
                      theta_star_coord2 = theta[2],
                      crossfit_cov = NA_real_)

# Create grid based on given theta values
results[, coord1 := rep(seq(theta[1] - 0.25, theta[1] + 0.25, by = 0.001), 
                        each = 501)]

results[, coord2 := rep(seq(theta[2] - 0.25, theta[2] + 0.25, by = 0.001), 
                        times = 501)]

# Set up progress bar
pb <- progress_bar$new(format = "sim :current / :total [:bar] :eta",
                       total = nrow(results), clear = T, show_after = 0)

# Generate n d-dimensional vectors Y. (Each observation is a row.)
set.seed(20200324)
Y <- mvrnorm(n = n_val, mu = theta, Sigma = diag(d_val))
 
# Set seed
set.seed(sim_val)

# Split Y into Y_0 and Y_1
Y_0_indices <- sample(1:n_val, size = n_val / 2)

Y_1_indices <- setdiff(1:n_val, Y_0_indices)

Y_0 <- Y[Y_0_indices, ]

Y_1 <- Y[Y_1_indices, ]

# Overall means of Y_0 and Y_1 observations
Y_0_bar <- apply(Y_0, MARGIN = 2, FUN = mean)

Y_1_bar <- apply(Y_1, MARGIN = 2, FUN = mean)

# Check whether each test point is within confidence region.
for(row in 1:nrow(results)) {
  
  # Increment progress bar
  pb$tick()
  
  # Store Y_0_bar and Y_1_bar vector
  results[row, Ybar0_coord1 := Y_0_bar[1]]
  results[row, Ybar0_coord2 := Y_0_bar[2]]
  results[row, Ybar1_coord1 := Y_1_bar[1]]
  results[row, Ybar1_coord2 := Y_1_bar[2]]
  
  # Extract coordinates for test point on grid
  test_point <- as.numeric(results[row, .(coord1, coord2)])
  
  # Does crossfit LRT universal confidence set cover true theta?
  results[row, crossfit_cov := 
            .5*(exp(-n_val * sum((Y_0_bar - Y_1_bar)^2) / 4 +
                      n_val * sum((Y_0_bar - test_point)^2) / 4) +
                  exp(-n_val * sum((Y_0_bar - Y_1_bar)^2) / 4 +
                        n_val * sum((Y_1_bar - test_point)^2) / 4)) <= 1/alpha]
  
}

# Save simulation results.
fwrite(results, 
       file = paste0("sim_data/fig01_crossfit_LRT_sim_",
                     sim_val, ".csv"))

