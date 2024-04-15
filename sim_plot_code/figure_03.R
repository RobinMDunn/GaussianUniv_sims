# Coverage/radius plots for split LRT around optimal proportions

library(tidyverse)
library(data.table)
library(grid)

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

# Read in data
radius_data <- fread(file = "data/radius/06_SLRT_opt_prop.csv")

# Get coverage proportions
coverage_data <- radius_data %>% 
  group_by(n, d, p0) %>% 
  summarise(coverage = mean(C_covered)) %>% 
  ungroup() %>% 
  mutate(dimension = factor(d, levels = c(10, 1000),
                            labels = c("d = 10", "d = 1000")))

sq_radius_data <- radius_data %>% 
  group_by(n, d, p0) %>% 
  summarise(mean_sq_rad = mean(C_sq_radius),
            sd_sq_rad = sd(C_sq_radius),
            n_sim = n()) %>% 
  ungroup() %>% 
  mutate(dimension = factor(d, levels = c(10, 1000),
                            labels = c("d = 10", "d = 1000")))

########################
##### Create plots #####
########################

# Optimal p0 plot
opt_p0 <- function(d, alpha) {
  return(1 + (2*d - sqrt(4*d^2 + 8*d*log(1/alpha))) / (4*log(1/alpha)))
}

d_vec <- seq(1, 250, by = 1)
p0_vec <- opt_p0(d_vec, alpha = 0.1)
p0_df <- data.frame(p0 = p0_vec, 
                    d = d_vec)

p0_plot <- p0_df %>% 
  ggplot(aes(x = d, y = p0)) +
  geom_line() +
  labs(x = "Dimension", 
       y = expression("Optimal"~p[0]),
       title = expression("Optimal"~p[0]~"by Dimension")) +
  paper_theme


# Coverage plot
coverage_plot <- 
  coverage_data %>% 
  ggplot(aes(x = p0, y = coverage)) +
  facet_wrap(dimension ~ ., scales = "free_x") +
  geom_point(alpha = 0.7) +
  geom_line(alpha = 0.5) +
  geom_hline(yintercept = 0.90, lty = "dashed") +
  labs(x = expression(p[0]), 
       y = "Coverage",
       title = expression("Coverage of Mult Normal Split LRT w/ Varying"~p[0])) +
  paper_theme

# Sq radius ratio plot
alpha <- 0.1

sq_radius_plot <- 
  sq_radius_data %>%
  mutate(expectation = 2*log(1/alpha)/(n*p0) + (1/(n*p0) + 1/(n*(1-p0)))*d,
         optimal = abs(opt_p0(d, alpha) - p0) < 0.01) %>% 
  ggplot(aes(x = p0, y = mean_sq_rad, col = optimal)) +
  facet_wrap(. ~ dimension, scales = "free") +
  geom_point(alpha = 0.8) +
  geom_errorbar(aes(ymin = (mean_sq_rad - qnorm(0.975)*sd_sq_rad) *
                      I(mean_sq_rad - qnorm(0.975)*sd_sq_rad > 0),
                    ymax = mean_sq_rad + qnorm(0.975)*sd_sq_rad),
                alpha = 1) +
  geom_line(aes(x = p0, y = expectation), alpha = 0.5, col = "red") +
  labs(x = expression(p[0]),
       y = "Squared radius",
       col = expression("Optimal"~p[0]),
       title = expression("Squared Radius of Multivariate Normal Split LRT w/ Varying"~p[0]),
       subtitle = "Red curve: expected squared radius") +
  scale_color_manual(values = c("blue", "darkgrey"), limits = c(T, F),
                     labels = c("Yes", "No")) +
  paper_theme

# Sq radius ratio plot, where color is not strictly necessary
sq_radius_plot_linetype <- sq_radius_data %>%
  mutate(expectation = 2*log(1/alpha)/(n*p0) + (1/(n*p0) + 1/(n*(1-p0)))*d,
         optimal = abs(opt_p0(d, alpha) - p0) < 0.01) %>% 
  ggplot(aes(x = p0, y = mean_sq_rad, col = optimal, linetype = optimal)) +
  facet_wrap(. ~ dimension, scales = "free") +
  geom_point(alpha = 0.8) +
  geom_errorbar(aes(ymin = (mean_sq_rad - qnorm(0.975)*sd_sq_rad) *
                      I(mean_sq_rad - qnorm(0.975)*sd_sq_rad > 0),
                    ymax = mean_sq_rad + qnorm(0.975)*sd_sq_rad),
                alpha = 1) +
  geom_line(aes(x = p0, y = expectation), alpha = 0.5, col = "red",
            inherit.aes = F) +
  labs(x = expression(p[0]),
       y = "Squared radius",
       col = expression("Optimal"~p[0]),
       linetype = expression("Optimal"~p[0]),
       title = expression("Squared Radius of Multivariate Normal Split LRT w/ Varying"~p[0]),
       subtitle = "Red curve: expected squared radius") +
  scale_linetype_manual(values = c("solid", "dashed"), limits = c(T, F),
                        labels = c("Yes", "No")) +
  scale_color_manual(values = c("blue", "darkgrey"), limits = c(T, F),
                     labels = c("Yes", "No")) +
  guides(color = guide_legend(override.aes = list(shape = NA))) +
  paper_theme

######################
##### Save plots #####
######################

ggsave(plot = p0_plot,
       filename = "plots/radius/06_p0_plot.pdf",
       width = 6, height = 4)

ggsave(plot = coverage_plot,
       filename = "plots/radius/06_coverage_plot.pdf",
       width = 8, height = 4)

ggsave(plot = sq_radius_plot,
       filename = "plots/radius/06_sq_radius_plot.pdf",
       width = 8, height = 4)

ggsave(plot = sq_radius_plot + labs(title = NULL, subtitle = NULL),
       filename = "plots/radius/06_sq_radius_plot_journal.pdf",
       width = 8, height = 3)

ggsave(plot = sq_radius_plot_linetype + labs(title = NULL, subtitle = NULL),
       filename = "plots/radius/06_sq_radius_plot_linetype_journal.pdf",
       width = 8, height = 3)
