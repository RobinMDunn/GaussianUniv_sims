# Simulations based on 
# "Comparing LRT to split SLR for the multivariate Normal case" section of 
# "Universal Inference Using the Split Likelihood Ratio Test."
# Save squared radius and coverage of A_n (usual confidence set)
# and C_n (universal confidence set).
# Exploring ratio of squared radii at small values of alpha.

suppressMessages(library(data.table))
suppressMessages(library(progress))
suppressMessages(library(MASS))

# Read in arguments for start/end d (dimension of Normal dist).
d_vec <- c(1, 2, 10, 100, 1000, 10000, 100000)
log_alpha_vec <- -10^(seq(0, 8, by = 0.5))
n_val <- 1000

# Construct data frame to store results
results <- data.table(n = n_val, 
                      expand.grid(log_alpha = log_alpha_vec, d = d_vec),
                      C_over_A_sq_rad_ub = NA_real_,
                      C_over_A_sq_rad_lb = NA_real_,
                      A_sq_radius = NA_real_,
                      C_sq_rad_expect = NA_real_,
                      true_ratio_expect = NA_real_)

# Set up progress bar
pb <- progress_bar$new(format = "sim :current / :total [:bar] :eta",
                       total = nrow(results), clear = T, show_after = 0)

# Simulate data, save sq radii of usual/split LRTs, and store ratio.
for(row in 1:nrow(results)) {
  
  # Increment progress bar
  pb$tick()
  
  # Extract n, d, alpha, and sim
  n_val <- results[row, n]
  
  d_val <- results[row, d]
  
  log_alpha_val <- results[row, log_alpha]
  
  # Get upper and lower bounds of ratio
  
  # Check condition from Prop 5.1 of Inglot (2010)
  if(d_val != 1 & log_alpha_val < log(0.17)) {
    
    # Upper bound on ratio of expectations
    results[row, 
            C_over_A_sq_rad_ub := (4*(-log_alpha_val)+ 4*d_val) / 
              (d_val + 2*(-log_alpha_val) - 5/2)]
    
  } 
  
  # Check condition derived from Inglot (2010) Thm A (due to Laurent/Massart) 
  # and Polland (2015) section 2.1 (due to Feller)
  if(d_val == 1 & log_alpha_val <= -10*(1+sqrt(5))/8) {
    
    # Upper bound on ratio of expectations
    results[row, 
            C_over_A_sq_rad_ub := (4*(-log_alpha_val)+ 4*d_val) / 
              (2*(-log_alpha_val) + 9 - 4*sqrt(5 - 2*log_alpha_val))]
    
  }
  
  # Lower bound on ratio of expectations. Derived from Inglot (2010) Thm A 
  # (due to Laurent/Massart).
  results[row, 
          C_over_A_sq_rad_lb := (4*(-log_alpha_val)+ 4*d_val) / 
            (2*(-log_alpha_val) + d_val + 2*sqrt(d*(-log_alpha_val)))]
  
  # Chi-sq cutoff
  c_alpha_d <- qchisq(p = log_alpha_val, df = d_val, lower.tail = F,
                      log.p = T)
  
  # Squared radius of A_n
  results[row, A_sq_radius := c_alpha_d / n_val]
  
  # Expected squared radius of C_n
  results[row, 
          C_sq_rad_expect := (4*(-log_alpha_val) + 4*d_val) / n_val]
  
  # Ratio of E[R^2(C_n) / R^2(A_n)]
  results[row, true_ratio_expect := C_sq_rad_expect / A_sq_radius]
}

fwrite(results, file = "sim_data/fig04.csv")
