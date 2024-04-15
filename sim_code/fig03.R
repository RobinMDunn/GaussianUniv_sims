# Simulations based on 
# "Comparing LRT to split SLR for the multivariate Normal case" section of 
# "Universal Inference Using the Split Likelihood Ratio Test."
# Save squared radius and coverage C_n (universal confidence set)
# varying proportions of observations in D_0 and D_1.
# Use range around optimized proportions.

suppressMessages(library(data.table))
suppressMessages(library(progress))
suppressMessages(library(MASS))

# Read in arguments for start/end d (dimension of Normal dist) and
# start/end sim indices.
start_d <- 10
end_d <- 1000
start_sim <- 1
end_sim <- 1000

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  args <- as.numeric(args)
  start_d <- args[1]
  end_d <- args[2]
  start_sim <- args[3]
  end_sim <- args[4]
}

# Set alpha level
alpha <- 0.1

# Construct vector of d values
d_all <- c(10, 1000)
d_vec <- d_all[start_d <= d_all & d_all <= end_d]

# Optimal p0 function
opt_p0 <- function(d, alpha) {
  return(1 + (2*d - sqrt(4*d^2 + 8*d*log(1/alpha))) / (4*log(1/alpha)))
}

# Construct vectors of p0 values
p0_d_mat <- matrix(NA, nrow = 0, ncol = 2)

if(10 %in% d_vec) {
  p0_d10 <- seq(from = opt_p0(d = 10, alpha = alpha) - 0.2, 
                to = opt_p0(d = 10, alpha = alpha) + 0.2,
                by = 0.025)
  
  p0_d_mat <- rbind(p0_d_mat, expand.grid(d = 10, p0 = p0_d10))
}

if(1000 %in% d_vec) {
  p0_d1000 <- seq(from = opt_p0(d = 1000, alpha = alpha) - 0.2, 
                 to = opt_p0(d = 1000, alpha = alpha) + 0.2,
                 by = 0.025)
  
  p0_d_mat <- rbind(p0_d_mat, expand.grid(d = 1000, p0 = p0_d1000))
}

# Number of simulations to perform at each combination of n and d
n_sim <- 1000

# Construct data frame to store results
results <- data.table(do.call(rbind, replicate(n_sim, p0_d_mat, simplify = F)),
                      n = 1000,
                      C_covered = NA_real_,
                      C_sq_radius = NA_real_) 

setorderv(results, cols = c("d", "p0"))

results[, sim := rep(1:n_sim, times = nrow(p0_d_mat))]

# Subset results for given simulations
results <- results[start_sim <= sim & sim <= end_sim]

# Set up progress bar
pb <- progress_bar$new(format = "sim :current / :total [:bar] :eta",
                       total = nrow(results), clear = T, show_after = 0)

# Simulate data, check whethersplit LRT at p0 contains true theta*, 
# and save sq radii.
for(row in 1:nrow(results)) {
  
  # Increment progress bar
  pb$tick()
  
  # Extract n, d, p0, and sim
  n_val <- results[row, n]
  
  d_val <- results[row, d]
  
  p0_val <- results[row, p0]
  
  sim_val <- results[row, sim]

  # Set theta vector
  theta <- rep(0, times = d_val)
  
  # Generate n d-dimensional vectors Y. (Each observation is a row.)
  Y <- t(apply(matrix(1:n_val, nrow = n_val), MARGIN = 1,
               FUN = function(x) rnorm(n = length(theta), 
                                       mean = theta, sd = 1)))
  
  # Split Y into Y_0 and Y_1
  Y_0_indices <- sample(1:n_val, size = round(n_val * p0_val))
  
  Y_1_indices <- setdiff(1:n_val, Y_0_indices)
  
  Y_0 <- Y[Y_0_indices, ]
  
  Y_1 <- Y[Y_1_indices, ]
  
  # Convert vectors to matrices with one row
  if(is.null(dim(Y_0))) {
    Y_0 <- matrix(Y_0, nrow = 1)
  }
  
  if(is.null(dim(Y_1))) {
    Y_1 <- matrix(Y_1, nrow = 1)
  }
  
  # Overall means of Y_0 and Y_1 observations
  Y_0_bar <- apply(Y_0, MARGIN = 2, FUN = mean)
  
  Y_1_bar <- apply(Y_1, MARGIN = 2, FUN = mean)
  
  # Does C_n (universal confidence set) cover true theta?
  results[row, C_covered := sum((theta - Y_0_bar)^2) <= 
            (2 / (n*p0_val)) * log(1/alpha) + sum((Y_0_bar - Y_1_bar)^2)]
  
  # Squared radius of C_n
  results[row, C_sq_radius := (2 / (n*p0_val)) * log(1/alpha) + 
            sum((Y_0_bar - Y_1_bar)^2)]
}

# Save simulation results
fwrite(results, 
       file = paste0("data/radius/06_SLRT_opt_prop_d_", start_d, "_", end_d, 
                     "_sim_", start_sim, "_", end_sim, ".csv"))