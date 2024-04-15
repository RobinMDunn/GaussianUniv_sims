# Empirical evaluation of r^2(C_n^split) / r^2(A_n) in multivariate normal case. 
# What proportion of times is this ratio <= 4?
# Or, equivalently, what proportion of the time is r^2(C_n^split) / r^2(A_n) <= 2?

library(data.table)
library(tidyverse)

# Read in data

radius_data_n10 <- fread(file = "data/radius/03_SLRT_sim10000_n_10.csv")

radius_data_n10_lowdim <- fread(file = "data/radius/03_SLRT_sim10000_n_10_d_1_20.csv")

radius_data_n1000 <- fread(file = "data/radius/03_SLRT_sim10000_n_1000.csv")

radius_data_n1000_lowdim <- fread(file = "data/radius/03_SLRT_sim10000_n_1000_d_1_20.csv")

alpha <- 0.1

# Create theme
paper_theme <- theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 14),
        legend.title = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        panel.spacing = unit(1.2, "lines"))

# Get proportion of simulations where squared radius ratio <= 4

prop_n10 <- radius_data_n10 %>% group_by(d) %>% 
  summarise(prop = mean(C_sq_radius / A_sq_radius <= 4)) %>% 
  mutate(upper_bound = 
           1 - alpha - 
           log(1/alpha) * dchisq(qchisq(p = 1 - alpha, df = d), df = d),
         lower_bound = 
           1 - alpha - 
           log(1/alpha) * dchisq(qchisq(p = 1 - alpha, df = d) + log(alpha),
                                 df = d),
         n = "10 observations")

prop_n10_lowdim <- radius_data_n10_lowdim %>% group_by(d) %>% 
  summarise(prop = mean(C_sq_radius / A_sq_radius <= 4)) %>% 
  mutate(upper_bound = 
           1 - alpha - 
           log(1/alpha) * dchisq(qchisq(p = 1 - alpha, df = d), df = d),
         lower_bound = 
           1 - alpha - 
           log(1/alpha) * dchisq(qchisq(p = 1 - alpha, df = d) + log(alpha),
                                 df = d),
         n = "10 observations")

prop_n1000 <- radius_data_n1000 %>% group_by(d) %>% 
  summarise(prop = mean(C_sq_radius / A_sq_radius <= 4)) %>% 
  mutate(upper_bound = 
           1 - alpha - 
           log(1/alpha) * dchisq(qchisq(p = 1 - alpha, df = d), df = d),
         lower_bound = 
           1 - alpha - 
           log(1/alpha) * dchisq(qchisq(p = 1 - alpha, df = d) + log(alpha),
                                 df = d),
         n = "1000 observations")

prop_n1000_lowdim <- radius_data_n1000_lowdim %>% group_by(d) %>% 
  summarise(prop = mean(C_sq_radius / A_sq_radius <= 4)) %>% 
  mutate(upper_bound = 
           1 - alpha - 
           log(1/alpha) * dchisq(qchisq(p = 1 - alpha, df = d), df = d),
         lower_bound = 
           1 - alpha - 
           log(1/alpha) * dchisq(qchisq(p = 1 - alpha, df = d) + log(alpha),
                                 df = d),
         n = "1000 observations")

all_prop <- rbind(prop_n10, prop_n1000)

all_prop_lowdim <- rbind(prop_n10_lowdim, prop_n1000_lowdim)

# Check that various computations of lower bound match
d <- 50
alpha <- 0.1
c_alpha_d <- qchisq(p = 1-alpha, df = d)

formula_calc <- log(1/alpha) * (1/2^(d/2)) * (1/gamma(d/2)) * 
  (c_alpha_d)^(d/2-1) * exp(-c_alpha_d/2)

dchisq_calc <- log(1/alpha) * dchisq(qchisq(p = 1 - alpha, df = d), df = d)

log_scale <- log(log(1/alpha)) - (d/2)*log(2) - lgamma(d/2) + 
  (d/2-1)*log(c_alpha_d) - (c_alpha_d/2)

exp_log_scale_calc <- exp(log_scale)

all.equal(formula_calc, dchisq_calc)

all.equal(formula_calc, exp_log_scale_calc)

# Get log(1/alpha)*f_d(c_{alpha, d}) and 
# log(1/alpha)*f_d(c_{alpha, d} + log(alpha)) at varying d and alpha
chisq_d_interp_df <- 
  expand.grid(alpha = c(0.01, 0.05, 0.10),
              d = seq(1, 1000, by = 1),
              lb_term = NA_real_,
              ub_term = NA_real_) %>% 
  as.data.table()

# Compute lower bound
chisq_d_interp_df[, lb_term := log(1/alpha) * 
                    dchisq(qchisq(p = 1 - alpha, df = d) + log(alpha), df = d),
                  by = seq_len(nrow(chisq_d_interp_df))]

# Compute upper bound
chisq_d_interp_df[, ub_term := log(1/alpha) * 
                    dchisq(qchisq(p = 1 - alpha, df = d), df = d),
                  by = seq_len(nrow(chisq_d_interp_df))]

# Plot log(1/alpha)*f_d(c_{alpha, d}) and 
# log(1/alpha)*f_d(c_{alpha, d} + log(alpha)) at varying d and alpha
chisq_d_interp_plot_dleq100 <- chisq_d_interp_df %>% 
  dplyr::mutate(alpha = factor(alpha, levels = c(0.1, 0.05, 0.01),
                               labels = c(0.1, 0.05, 0.01))) %>% 
  dplyr::filter(d <= 100) %>% 
  pivot_longer(cols = c(lb_term, ub_term)) %>% 
  dplyr::mutate(name = factor(name, levels = c("lb_term", "ub_term"),
                              labels = c("log(1/alpha) ~ f[d](c[alpha][','][d] + log(alpha))",
                                         "log(1/alpha) ~ f[d](c[alpha][','][d])"))) %>% 
  ggplot(aes(x = d, y = value, col = as.factor(alpha))) +
  facet_wrap(. ~ name, labeller = label_parsed) +
  geom_line() +
  paper_theme +
  labs(col = expression(alpha),
       x = "Dimension",
       y = "Value") +
  scale_color_manual(values = c("blue", "red", "black"))

chisq_d_interp_plot_dleq20 <- chisq_d_interp_df %>% 
  dplyr::mutate(alpha = factor(alpha, levels = c(0.1, 0.05, 0.01),
                               labels = c(0.1, 0.05, 0.01))) %>% 
  dplyr::filter(d <= 20) %>% 
  pivot_longer(cols = c(lb_term, ub_term)) %>% 
  dplyr::mutate(name = factor(name, levels = c("lb_term", "ub_term"),
                              labels = c("log(1/alpha) ~ f[d](c[alpha][','][d] + log(alpha))",
                                         "log(1/alpha) ~ f[d](c[alpha][','][d])"))) %>% 
  ggplot(aes(x = d, y = value, col = as.factor(alpha))) +
  facet_wrap(. ~ name, labeller = label_parsed) +
  geom_line() +
  paper_theme +
  labs(col = expression(alpha),
       x = "Dimension",
       y = "Value") +
  scale_color_manual(values = c("blue", "red", "black"))
  
chisq_d_interp_plot_dleq1000 <- chisq_d_interp_df %>% 
  dplyr::mutate(alpha = factor(alpha, levels = c(0.1, 0.05, 0.01),
                               labels = c(0.1, 0.05, 0.01))) %>% 
  dplyr::filter(d <= 1000) %>% 
  pivot_longer(cols = c(lb_term, ub_term)) %>% 
  dplyr::mutate(name = factor(name, levels = c("lb_term", "ub_term"),
                              labels = c("log(1/alpha) ~ f[d](c[alpha][','][d] + log(alpha))",
                                         "log(1/alpha) ~ f[d](c[alpha][','][d])"))) %>% 
  ggplot(aes(x = d, y = value, col = as.factor(alpha))) +
  facet_wrap(. ~ name, labeller = label_parsed) +
  geom_line() +
  paper_theme +
  labs(col = expression(alpha),
       x = "Dimension",
       y = "Value") +
  scale_color_manual(values = c("blue", "red", "black"))

#################
##### Plots #####
#################

prop_plot_highdim <- all_prop %>% 
  pivot_longer(cols = c( "upper_bound", "lower_bound"), 
               names_to = "bound_label", values_to = "bound_value") %>% 
  mutate(bound_label = factor(bound_label, 
                              levels = c("upper_bound", "lower_bound"),
                              labels = c("Upper bound", "Lower bound"))) %>% 
  ggplot(aes(x = d, y = bound_value, color = bound_label)) +
  facet_grid(. ~ n) +
  geom_point(aes(x = d, y = prop), col = "black", alpha = 0.3) +
  geom_line(size = 1, alpha = 0.7) +
  scale_color_manual(values = c("red", "blue")) +
  labs(x = "Dimension",
       y = "Proportion",
       color = expression(atop("Bound on", P(r^{2}~(C[n]^split*(alpha))~
                            "/"~r^{2}~(C[n]^LRT*(alpha))<="4"))),
       title = expression(atop("Proportion of Simulations where"~r^{2}~(C[n]^split*(alpha))~
                                 "/"~r^{2}~(C[n]^LRT*(alpha))<="4",
                               "With Lower and Upper Bounds on Probability")),
       subtitle = "Higher dimensions (d = 10 to d = 1000)") +
  paper_theme +
  theme(legend.title = element_text(size = 12))

prop_plot_lowdim <- all_prop_lowdim %>% 
  pivot_longer(cols = c( "upper_bound", "lower_bound"), 
               names_to = "bound_label", values_to = "bound_value") %>% 
  mutate(bound_label = factor(bound_label, 
                              levels = c("upper_bound", "lower_bound"),
                              labels = c("Upper bound", "Lower bound"))) %>% 
  ggplot(aes(x = d, y = bound_value, color = bound_label)) +
  facet_grid(. ~ n) +
  geom_point(aes(x = d, y = prop), col = "black", alpha = 0.3) +
  geom_line(size = 1, alpha = 0.7) +
  scale_color_manual(values = c("red", "blue")) +
  labs(x = "Dimension",
       y = "Proportion",
       color = expression(atop("Bound on", P(r^{2}~(C[n]^split*(alpha))~
                                               "/"~r^{2}~(C[n]^LRT*(alpha))<="4"))),
       title = expression(atop("Proportion of Simulations where"~r^{2}~(C[n]^split*(alpha))~
                            "/"~r^{2}~(C[n]^LRT*(alpha))<="4",
                            "With Lower and Upper Bounds on Probability")),
       subtitle = "Lower dimensions (d = 1 to d = 20)") +
  paper_theme +
  theme(legend.title = element_text(size = 12))

prop_plot_d1to20 <- all_prop_lowdim %>% 
  pivot_longer(cols = c( "upper_bound", "lower_bound"), 
               names_to = "bound_label", values_to = "bound_value") %>% 
  mutate(bound_label = factor(bound_label, 
                              levels = c("upper_bound", "lower_bound"),
                              labels = c("Upper bound", "Lower bound"))) %>% 
  ggplot(aes(x = d, y = bound_value, color = bound_label)) +
  facet_grid(. ~ n) +
  geom_line(size = 1, alpha = 0.7) +
  geom_point(aes(x = d, y = prop), col = "black", alpha = 0.3) +
  scale_color_manual(values = c("red", "blue")) +
  scale_y_continuous(limits = c(-0.3, 1),
                     breaks = c(-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1)) +
  geom_hline(yintercept = 0.9, lty = "dashed") +
  labs(x = "Dimension",
       y = "Proportion",
       color = expression(atop("Bound on", P(r(C[n]^split*(alpha))~
                                               "/"~r(C[n]^LRT*(alpha))<="2"))),
       title = expression(atop("Proportion of Simulations where"~r*
                                 (C[n]^split*(alpha))~
                                 "/"~r(C[n]^LRT*(alpha))<="2",
                               "With Lower and Upper Bounds on Probability"))) +
  paper_theme +
  theme(legend.title = element_text(size = 12))

prop_plot_d1to20_linetype <- all_prop_lowdim %>% 
  pivot_longer(cols = c( "upper_bound", "lower_bound"), 
               names_to = "bound_label", values_to = "bound_value") %>% 
  mutate(bound_label = factor(bound_label, 
                              levels = c("upper_bound", "lower_bound"),
                              labels = c("Upper bound", "Lower bound"))) %>% 
  ggplot(aes(x = d, y = bound_value, color = bound_label,
             linetype = bound_label)) +
  facet_grid(. ~ n) +
  geom_line(size = 1, alpha = 0.7) +
  geom_point(aes(x = d, y = prop), col = "black", alpha = 0.3) +
  scale_color_manual(values = c("red", "blue")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_y_continuous(limits = c(-0.3, 1),
                     breaks = c(-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1)) +
  geom_hline(yintercept = 0.9, lty = "dashed") +
  labs(x = "Dimension",
       y = "Proportion",
       color = expression(atop("Bound on", P(r(C[n]^split*(alpha))~
                                               "/"~r(C[n]^LRT*(alpha))<="2"))),
       linetype = expression(atop("Bound on", P(r(C[n]^split*(alpha))~
                                               "/"~r(C[n]^LRT*(alpha))<="2"))),
       title = expression(atop("Proportion of Simulations where"~r*
                                 (C[n]^split*(alpha))~
                                 "/"~r(C[n]^LRT*(alpha))<="2",
                               "With Lower and Upper Bounds on Probability"))) +
  paper_theme +
  theme(legend.title = element_text(size = 12),
        legend.key.width = unit(2, "cm"))

######################
##### Save plots #####
######################

ggsave(plot = prop_plot_highdim,
       filename = "plots/radius/10_ratio_leq_4_highdim.pdf",
       width = 8.2, height = 4)

ggsave(plot = prop_plot_lowdim,
       filename = "plots/radius/10_ratio_leq_4_lowdim.pdf",
       width = 8.2, height = 4)

ggsave(plot = prop_plot_d1to20,
       filename = "plots/radius/10_ratio_leq_4_d1to20.pdf",
       width = 8.2, height = 4)

ggsave(plot = prop_plot_d1to20 + labs(title = NULL, subtitle = NULL),
       filename = "plots/radius/10_ratio_leq_4_d1to20_journal.pdf",
       width = 8.2, height = 3)

ggsave(plot = prop_plot_d1to20_linetype + labs(title = NULL, subtitle = NULL),
       filename = "plots/radius/10_ratio_leq_4_d1to20_linetype_journal.pdf",
       width = 8.2, height = 3)

ggsave(plot = chisq_d_interp_plot_dleq100 + labs(title = NULL, subtitle = NULL),
       filename = "plots/radius/10_chisq_d_interp_plot_dleq100.pdf",
       width = 8.2, height = 4)

ggsave(plot = chisq_d_interp_plot_dleq20 + labs(title = NULL, subtitle = NULL),
       filename = "plots/radius/10_chisq_d_interp_plot_dleq20.pdf",
       width = 8.2, height = 3)

ggsave(plot = chisq_d_interp_plot_dleq1000 + labs(title = NULL, subtitle = NULL),
       filename = "plots/radius/10_chisq_d_interp_plot_dleq1000.pdf",
       width = 8.2, height = 4)
