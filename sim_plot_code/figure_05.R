# Empirical evaluation of r^2(C_n^split) / r^2(A_n) in multivariate normal case. 
# What proportion of times is this ratio <= 4?
# Or, equivalently, what proportion of the time is r^2(C_n^split) / r^2(A_n) <= 2?

library(data.table)
library(tidyverse)

# Read in data
radius_data_n10_lowdim <- fread(file = "sim_data/fig05_n_10_d_1_20_by_1.csv")

radius_data_n1000_lowdim <- fread(file = "sim_data/fig05_n_1000_d_1_20_by_1.csv")

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

all_prop_lowdim <- rbind(prop_n10_lowdim, prop_n1000_lowdim)

#################
##### Plots #####
#################

prop_plot_lowdim <- all_prop_lowdim %>% 
  pivot_longer(cols = c( "upper_bound", "lower_bound"), 
               names_to = "bound_label", values_to = "bound_value") %>% 
  mutate(bound_label = factor(bound_label, 
                              levels = c("upper_bound", "lower_bound"),
                              labels = c("Upper bound", "Lower bound"))) %>% 
  ggplot(aes(x = d, y = bound_value, color = bound_label)) +
  facet_grid(. ~ n) +
  geom_point(aes(x = d, y = prop), col = "black", alpha = 0.3) +
  geom_line(linewidth = 1, alpha = 0.7) +
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

######################
##### Save plots #####
######################

ggsave(plot = prop_plot_lowdim + labs(title = NULL, subtitle = NULL),
       filename = "sim_plots/figure_05.pdf",
       width = 8.2, height = 3)
