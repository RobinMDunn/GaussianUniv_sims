# Simulations based on 
# "Comparing LRT to split SLR for the multivariate Normal case" section of 
# "Universal Inference Using the Split Likelihood Ratio Test."
# Save squared radius and coverage of A_n (usual confidence set)
# and C_n (universal confidence set).

suppressMessages(library(data.table))
suppressMessages(library(progress))

# Read in arguments for start/end d (dimensional of Normal dist).
start_d <- 10
end_d <- 1000
n_input <- 1000
increment <- 10

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  args <- as.numeric(args)
  start_d <- args[1]
  end_d <- args[2]
  n_input <- args[3]
  increment <- args[4]
}

# Set alpha level
alpha <- 0.1

# Construct vector of d values
d_vec <- seq(start_d, end_d, by = increment)

# Number of simulations to perform at each combination of n and d
n_sim <- 10000

# Construct data frame to store results
results <- data.table(n = n_input, 
                      expand.grid(d = d_vec, sim = 1:n_sim),
                      A_sq_radius = NA_real_,
                      C_sq_radius = NA_real_)

# Set up progress bar
pb <- progress_bar$new(format = "sim :current / :total [:bar] :eta",
                       total = nrow(results), clear = T, show_after = 0)

# Simulate data and save sq radii of usual/split LRT.
for(row in 1:nrow(results)) {
  
  # Increment progress bar
  pb$tick()
  
  # Extract n, d, and sim
  n_val <- results[row, n]
  
  d_val <- results[row, d]
  
  sim_val <- results[row, sim]
  
  # Set seed
  set.seed(7*n_val + 8*d_val + sim_val)
  
  # Simulate theta vector
  theta <- runif(n = d_val, min = 0, max = 1)
  
  # Generate n d-dimensional vectors Y. (Each observation is a row.)
  Y <- t(apply(matrix(1:n_val, nrow = n_val), MARGIN = 1,
                FUN = function(x) rnorm(n = length(theta), 
                                        mean = theta, sd = 1)))

  
  # Overall mean of Y observations
  Y_bar <- apply(Y, MARGIN = 2, FUN = mean)
  
  # Chi-sq cutoff
  # c_alpha_d <- qchisq(p = 1 - alpha, df = d_val, lower.tail = T)
  
  # Does A_n (usual confidence set) cover true theta?
  # results[row, A_covered := sum((theta - Y_bar)^2) <= c_alpha_d / n_val]
  
  # Squared radius of A_n
  results[row, A_sq_radius := qchisq(p = 1 - alpha, df = d_val,
                                     lower.tail = T) / n_val]
  
  # Split Y into Y_0 and Y_1
  Y_0_indices <- sample(1:n_val, size = n_val / 2)
  
  Y_1_indices <- setdiff(1:n_val, Y_0_indices)
  
  if(d_val > 1) {
    
    Y_0 <- Y[Y_0_indices, ]
    
    Y_1 <- Y[Y_1_indices, ]
    
    # Overall means of Y_0 and Y_1 observations
    Y_0_bar <- apply(Y_0, MARGIN = 2, FUN = mean)
    
    Y_1_bar <- apply(Y_1, MARGIN = 2, FUN = mean)
    
  } else {
    
    Y_0 <- Y[Y_0_indices]
    
    Y_1 <- Y[Y_1_indices]
    
    # Overall means of Y_0 and Y_1 observations
    Y_0_bar <- mean(Y_0)
    
    Y_1_bar <- mean(Y_1)
    
  }
  
  
  # Does C_n (universal confidence set) cover true theta?
  # results[row, C_covered := sum((theta - Y_0_bar)^2) <= 
  #           4*log(1/alpha) / n_val + sum((Y_0_bar - Y_1_bar)^2)]
  
  # Squared radius of C_n
  results[row, 
          C_sq_radius := 4*log(1/alpha) / n_val + sum((Y_0_bar - Y_1_bar)^2)]
}

# Save simulation results.
fwrite(results, 
       file = paste0("data/radius/03_SLRT_sim10000_n_", n_input, "_d_", 
                     start_d, "_", end_d, ".csv"))