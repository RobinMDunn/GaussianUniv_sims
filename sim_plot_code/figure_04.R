# E[r^2(C_n^split) / r^2(A_n)] at small alpha for varying d 

library(data.table)
library(tidyverse)

# Read in data

small_alpha_data <- fread(file = "data/radius/04_SLRT_small_alpha.csv")

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

# Plot expected ratio of squared radii (and bounds) for six d values
ratio_plot <- small_alpha_data %>% 
  filter(d <= 10000, log_alpha <= -10) %>% 
  pivot_longer(cols = c(C_over_A_sq_rad_lb, C_over_A_sq_rad_ub,
                        true_ratio_expect),
               names_to = "ratio_label",
               values_to = "ratio_value") %>% 
  mutate(d = factor(d, levels = c(1, 2, 10, 100, 1000, 10000),
                    labels = c("d = 1", "d = 2", "d = 10", "d = 100",
                               "d = 1,000", "d = 10,000")),
         ratio_label = factor(ratio_label, 
                              levels = c("C_over_A_sq_rad_ub",
                                         "true_ratio_expect",
                                         "C_over_A_sq_rad_lb"),
                              labels = c("Upper bound", "True value",
                                         "Lower bound"))) %>% 
  ggplot(aes(x = log(-log_alpha, base = 10), y = ratio_value, col = ratio_label)) +
  facet_wrap(. ~ d, nrow = 2) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.7) +
  labs(x = expression(alpha),
       y = "Ratio of squared radii",
       title = expression("Bounds and Value of E["*r^{2}~(C[n]^{split}*(alpha))*
                            "] /"~r^{2}~(A[n](alpha))),
       col = "Ratio") +
  scale_x_continuous(breaks = c(2, 4, 6, 8),
                     labels = c(expression("e"**{-10**{2}}), 
                                expression("e"**{-10**{4}}),
                                expression("e"**{-10**{6}}),
                                expression("e"**{-10**{8}})),
                     trans = "reverse") +
  scale_color_manual(values = c("red", "black", "blue")) +
  paper_theme +
  theme(legend.position = "bottom")

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
       title = expression("Bounds and Value of E["*r^{2}~(C[n]^{split}*(alpha))*
                            "] /"~r^{2}~(C[n]^{LRT}*(alpha))),
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

# Plot expected ratio of squared radii (and bounds) for d = 10 and d = 100,000.
# Use different line types.
ratio_small_alpha_two_panel_linetype <- small_alpha_data %>% 
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
  ggplot(aes(x = log(-log_alpha, base = 10), y = ratio_value, col = ratio_label,
             linetype = ratio_label)) +
  facet_wrap(. ~ d, nrow = 1) +
  geom_point(alpha = 0.5) +
  geom_line(alpha = 0.7) +
  labs(x = expression(alpha),
       y = "Ratio of squared radii",
       title = expression("Bounds and Value of E["*r^{2}~(C[n]^{split}*(alpha))*
                            "] /"~r^{2}~(C[n]^{LRT}*(alpha))),
       col = "Ratio",
       linetype = "Ratio") +
  scale_x_continuous(breaks = c(0, 2, 4, 6, 8),
                     labels = c(expression("e"**{-10**{0}}),
                                expression("e"**{-10**{2}}), 
                                expression("e"**{-10**{4}}),
                                expression("e"**{-10**{6}}),
                                expression("e"**{-10**{8}})),
                     trans = "reverse") +
  scale_color_manual(values = c("red", "black", "blue")) +
  scale_linetype_manual(values = c("longdash", "solid", "dashed")) +
  paper_theme +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(shape = NA)))

######################
##### Save plots #####
######################

ggsave(plot = ratio_plot,
       filename = "plots/radius/04_small_alpha.pdf",
       width = 7, height = 5.5)

ggsave(plot = ratio_small_alpha_two_panel,
       filename = "plots/radius/04_small_alpha_2panel.pdf",
       width = 7, height = 3.5)

ggsave(plot = ratio_small_alpha_two_panel + labs(title = NULL, subtitle = NULL),
       filename = "plots/radius/04_small_alpha_2panel_journal.pdf",
       width = 7, height = 3)

ggsave(plot = ratio_small_alpha_two_panel_linetype + labs(title = NULL, subtitle = NULL),
       filename = "plots/radius/04_small_alpha_2panel_linetype_journal.pdf",
       width = 7, height = 3)

