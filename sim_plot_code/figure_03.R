# Squared radius plot for split LRT around optimal proportions

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
radius_data <- fread(file = "sim_data/fig03.csv")

# Get squared radius summaries
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

# Sq radius ratio plot
alpha <- 0.1

sq_radius_plot <- sq_radius_data %>%
  mutate(expectation = 2*log(1/alpha)/(n*p0) + (1/(n*p0) + 1/(n*(1-p0)))*d,
         optimal = abs(opt_p0(d, alpha) - p0) < 0.01,
         ymin_val = as.numeric((mean_sq_rad - qnorm(0.975)*sd_sq_rad) *
                                 I(mean_sq_rad - qnorm(0.975)*sd_sq_rad > 0)),
         ymax_val = as.numeric(mean_sq_rad + qnorm(0.975)*sd_sq_rad)) %>% 
  ggplot(aes(x = p0, y = mean_sq_rad, col = optimal)) +
  geom_point(alpha = 0.8) +
  facet_wrap(. ~ dimension, scales = "free") +
  geom_errorbar(aes(x = p0, ymin = ymin_val, ymax = ymax_val),
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

#####################
##### Save plot #####
#####################

ggsave(plot = sq_radius_plot + labs(title = NULL, subtitle = NULL),
       filename = "sim_plots/figure_03.pdf",
       width = 8, height = 3)
