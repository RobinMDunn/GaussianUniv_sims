# E[r^2(C_n^split) / r^2(A_n)] at small alpha for varying d 

library(data.table)
library(tidyverse)

# Read in data

small_alpha_data <- fread(file = "sim_data/fig04.csv")

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

#################
##### Plots #####
#################

# Plot expected ratio of squared radii (and bounds) for d = 10 and d = 100,000
ratio_small_alpha_two_panel <- small_alpha_data %>% 
  filter(d %in% c(10, 100000)) %>% 
  pivot_longer(cols = c(C_over_A_sq_rad_lb, C_over_A_sq_rad_ub,
                        true_ratio_expect),
               names_to = "ratio_label",
               values_to = "ratio_value") %>% 
  mutate(d = factor(d, levels = c(10, 100000),
                    labels = c("d = 10", "d = 100,000")),
         ratio_label = factor(ratio_label, 
                              levels = c("C_over_A_sq_rad_ub",
                                         "true_ratio_expect",
                                         "C_over_A_sq_rad_lb"),
                              labels = c("Upper bound", "True value",
                                         "Lower bound"))) %>% 
  ggplot(aes(x = log(-log_alpha, base = 10), y = ratio_value, col = ratio_label)) +
  facet_wrap(. ~ d, nrow = 1) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.7) +
  labs(x = expression(alpha),
       y = "Ratio of squared radii",
       col = "Ratio") +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8),
                     labels = c(expression("e"**{-10**{0}}),
                                expression("e"**{-10**{2}}), 
                                expression("e"**{-10**{4}}),
                                expression("e"**{-10**{6}}),
                                expression("e"**{-10**{8}})),
                     trans = "reverse") +
  scale_color_manual(values = c("red", "black", "blue")) +
  paper_theme +
  theme(legend.position = "bottom")

######################
##### Save plots #####
######################

ggsave(plot = ratio_small_alpha_two_panel,
       filename = "sim_plots/figure_04.pdf",
       width = 7, height = 3.5)
