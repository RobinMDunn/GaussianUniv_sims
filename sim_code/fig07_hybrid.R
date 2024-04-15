# Simulate power of subsampling hybrid LRT.
# Within subsamples, test statistics are 
# R_n I(||\bar{Y}_1|| > 1) + T_n(\hat{\theta}_0) I(||\bar{Y}_1|| < 0.5) +
#   1 I(||\bar{Y}_1|| \in [0.5, 1]).
# Use split statistic of
# T_n(\hat{\theta}_0) = L_0(\hat{\theta}_1) / L_0(\hat{\theta}_0) with choice
# \hat{\theta}_1 same as RIPR numerator.
#
# H_0: ||theta*|| \in [0.5, 1.0]
# H_1: ||theta*|| \notin [0.5, 1.0]

suppressMessages(library(data.table))
suppressMessages(library(progress))
suppressMessages(library(MASS))
suppressMessages(library(tidyverse))

# Read in arguments for d (dimensional of Normal dist),
# start coordinate, end coordinate, increment, 
# number of simulations, number of subsamples, number of observations,
# and whether to compute the test stat (0 is faster, but 1 computes test stat).

d_val <- 2
coord_start <- 0
coord_end <- 2
increment <- 0.05
n_sim <- 1000
B <- 100
n_obs <- 1000
compute_ts <- 0

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  args <- as.numeric(args)
  d_val <- args[1]
  coord_start <- args[2]
  coord_end <- args[3]
  increment <- args[4]
  n_sim <- args[5]
  B <- args[6]
  n_obs <- args[7]
  compute_ts <- args[8]
}

# Set alpha level
alpha_val <- 0.1

# Construct data frame to store results
results <- data.table(n = n_obs, 
                      n_sim = n_sim,
                      B = B,
                      alpha = alpha_val,
                      d = d_val,
                      expand.grid(sim = 1:n_sim,
                                  theta_star_1 = seq(coord_start, coord_end, 
                                                     by = increment)),
                      test_stat = NA_real_,
                      reject = NA_real_,
                      case_1_prop = NA_real_, 
                      case_2_prop = NA_real_,
                      case_3_prop = NA_real_)

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
  
  # Vector for subsampled likelihood ratios
  lik_ratio <- rep(NA, B)
  
  # Variables for case 1/2/3 proportions
  case_1 <- case_2 <- case_3 <- 0
  
  # Generate likelihood ratios through subsampling
  for(subsamp in 1:B) {
    
    # Split Y into Y_0 and Y_1
    Y_0_indices <- sample(1:n_obs, size = n_obs / 2)
    
    Y_1_indices <- setdiff(1:n_obs, Y_0_indices)
    
    Y_0 <- matrix(Y[Y_0_indices, ], ncol = d_val)
    
    Y_1 <- matrix(Y[Y_1_indices, ], ncol = d_val)
    
    # Determine ||\bar{Y}_1|| to choose method
    Y_bar_1 <- apply(Y_1, MARGIN = 2, FUN = mean)
    
    norm_Y_bar_1 <- sqrt(sum(Y_bar_1^2))
     
    # Compute test stat based on value of ||\bar{Y}_1||
    if(norm_Y_bar_1 >= 0.5 & norm_Y_bar_1 <= 1) {
      
      # ||\bar{Y}_1|| satisfies case 2
      case_2 <- case_2 + 1/B
      
      # If ||\bar{Y}_1|| \in [0.5, 1], set lklihood ratio for the subsample to 1
      lik_ratio[subsamp] <- 1
      
    } else if(norm_Y_bar_1 > 1) {
      
      # ||\bar{Y}_1|| satisfies case 3
      case_3 <- case_3 + 1/B
      
      # Use RIPR stat if ||\bar{Y}_1|| > 1
      
      # Set theta_hat_1 to \bar{Y}_1.
      theta_hat_1 <- Y_bar_1
      
      # Set theta_hat_0 as RIPR projection of p_{theta_hat_1} onto 
      # p_theta for theta\in\Theta_0.
      # theta_hat_0 = theta_hat_1 / ||theta_hat_1|| if ||theta_hat_1|| > 1
      theta_hat_0 <- theta_hat_1 / norm_Y_bar_1
      
      # Compute likelihood
      lik_ratio[subsamp] <- exp(
        -0.5 * sum(apply(Y_0, MARGIN = 1, FUN = function(row) sum((row - theta_hat_1)^2))) + 
          0.5 * sum(apply(Y_0, MARGIN = 1, FUN = function(row) sum((row - theta_hat_0)^2)))
      )
      
    } else if(norm_Y_bar_1 < 0.5) {
      
      # ||\bar{Y}_1|| satisfies case 1
      case_1 <- case_1 + 1/B
      
      # Use split LRT if ||\bar{Y}_1|| < 0.5
      
      # Set theta_hat_1 to \bar{Y}_1.
      theta_hat_1 <- Y_bar_1
      
      # Set theta_hat_0 as MLE under H_0, constructed from D_0.
      # theta_hat_0 = \bar{Y}_0 if ||\bar{Y}_0|| \in [0.5, 1.0]
      # theta_hat_0 = \bar{Y}_0 / ||\bar{Y}_0|| if ||\bar{Y}_0|| > 1.0
      # theta_hat_0 = 0.5* \bar{Y}_0 / ||\bar{Y}_0|| if ||\bar{Y}_0|| < 0.5
      theta_hat_0 <- apply(Y_0, MARGIN = 2, FUN = mean)
      
      theta_hat_0_norm <- sqrt(sum(theta_hat_0^2))
      
      if(theta_hat_0_norm > 1) {
        theta_hat_0 <- theta_hat_0 / theta_hat_0_norm
      } else if(theta_hat_0_norm < 0.5) {
        theta_hat_0 <- 0.5 * theta_hat_0 / theta_hat_0_norm
      }
      
      # Compute likelihood
      lik_ratio[subsamp] <- exp(
        -0.5 * sum(apply(Y_0, MARGIN = 1, FUN = function(row) sum((row - theta_hat_1)^2))) + 
          0.5 * sum(apply(Y_0, MARGIN = 1, FUN = function(row) sum((row - theta_hat_0)^2)))
      )
     
    }
    
    # Stop early if likelihood ratio already guarantees rejection
    if(compute_ts == 0 & 
       sum(lik_ratio[1:subsamp]) > B / alpha_val) {
      lik_ratio[(subsamp + 1):B] <- 0
      break
    }
  
  }
  
  # Store test statistic and case 1/2/3 proportions
  if(compute_ts == 1) {
    results[row, test_stat := mean(lik_ratio)]
    
    results[row, case_1_prop := case_1]
    results[row, case_2_prop := case_2]
    results[row, case_3_prop := case_3]
  }
  
  # Store decision. 
  results[row, reject := as.numeric(mean(lik_ratio) >= 1 / alpha_val)]

}

# Save rejection proportions
results <- results %>% 
  group_by(n, n_sim, B, alpha, d, theta_star_1) %>% 
  dplyr::summarise(reject = mean(reject),
                   mean_test_stat = mean(test_stat),
                   mean_case_1_prop = mean(case_1_prop),
                   mean_case_2_prop = mean(case_2_prop),
                   mean_case_3_prop = mean(case_3_prop))

# Save simulation results.
fwrite(results, file = paste0("data/RIPR/norm_interval/03_hybrid_d_",
                              d_val, "_n_", n_obs, "_coord_", 
                              coord_start, ".csv"))

