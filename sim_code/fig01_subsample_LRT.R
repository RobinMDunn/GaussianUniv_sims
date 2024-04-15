# Simulations based on 
# "Comparing LRT to split LRT for the multivariate Normal case" section of 
# "Universal Inference Using the Split Likelihood Ratio Test."
# Create coverage grid for subsampling split LRT.

suppressMessages(library(data.table))
suppressMessages(library(progress))
suppressMessages(library(MASS))
suppressMessages(library(Rcpp))

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

# Number of times to subsample
B <- 100

# Function to compute test statistic
suppressMessages(cppFunction('double test_stat_cpp(double n_val, 
                                                   NumericVector theta, 
                                                   NumericMatrix Y_0_bar, 
                                                   NumericMatrix Y_1_bar) {
  double sum_test_stat = 0;
  
  int B = Y_0_bar.nrow();
  
  int d_val = Y_0_bar.ncol();
  
  double sq_norm_Y0_Y1_bar;
  
  double sq_norm_Y0_bar_theta;
  
  for(int b = 0; b < B; b++) {
  
    sq_norm_Y0_Y1_bar = 0;
    
    sq_norm_Y0_bar_theta = 0;
    
    for(int dim = 0; dim < d_val; dim++) {
    
      sq_norm_Y0_Y1_bar += (Y_0_bar(b, dim) - Y_1_bar(b, dim)) * 
                             (Y_0_bar(b, dim) - Y_1_bar(b, dim));
    
      sq_norm_Y0_bar_theta += (Y_0_bar(b, dim) - theta[dim]) *
                                (Y_0_bar(b, dim) - theta[dim]);
      
    }
    
    sum_test_stat += exp((-n_val/4) * sq_norm_Y0_Y1_bar + 
                            (n_val/4) * sq_norm_Y0_bar_theta);
                            
  }
  
  double mean_test_stat = sum_test_stat / double(B);
  
  return mean_test_stat;
}'))

# Construct data frame to store results
results <- data.table(n = n_val, 
                      d = d_val,
                      coord1 = NA_real_,
                      coord2 = NA_real_,
                      sim = rep(sim_val, each = 501 * 501),
                      theta_star_coord1 = theta[1],
                      theta_star_coord2 = theta[2],
                      covered = NA_real_)

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

# Create matrices to store Y_0_bar and Y_1_bar
Y_0_bar <- matrix(NA, nrow = B, ncol = d_val)
Y_1_bar <- matrix(NA, nrow = B, ncol = d_val)

# Get B subsamples
for(subsample in 1:B) {
  
  # Split Y into Y_0 and Y_1
  Y_0_indices <- sample(1:n_val, size = n_val / 2)
  
  Y_1_indices <- setdiff(1:n_val, Y_0_indices)
  
  Y_0 <- as.matrix(Y[Y_0_indices, ])
  
  Y_1 <- as.matrix(Y[Y_1_indices, ])
  
  # Overall means of Y_0 and Y_1 observations
  Y_0_bar[subsample, ] <- apply(Y_0, MARGIN = 2, FUN = mean)
  
  Y_1_bar[subsample, ] <- apply(Y_1, MARGIN = 2, FUN = mean)
  
}

# Check whether each test point is within confidence region.
for(row in 1:nrow(results)) {
  
  # Increment progress bar
  pb$tick()
  
  # Extract coordinates for test point on grid
  test_point <- as.numeric(results[row, .(coord1, coord2)])
  
  # Does C_n (universal confidence set) cover true theta?
  results[row, 
          covered :=  test_stat_cpp(n_val = n_val, theta = test_point, 
                                    Y_0_bar = Y_0_bar, 
                                    Y_1_bar = Y_1_bar) <= 1 / alpha]
}

# Save simulation results.
fwrite(results, 
       file = paste0("data/region/fixed_theta_fixed_data/",
                     "03_normal_subsample_sim_",
                     sim_val, ".csv"))

