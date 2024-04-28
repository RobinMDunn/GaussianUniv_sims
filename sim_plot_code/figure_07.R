# Plot power and validity of subsampled split, subsampled hybrid RIPR/split, and 
# intersection tests.
# H_0: ||theta*|| \in [0.5, 1.0]
# H_1: ||theta*|| \notin [0.5, 1.0]

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

results_intersect <- fread("sim_data/fig07_intersect.csv")
results_intersect[, Method := "1. Intersection"]

results_split <- fread("sim_data/fig07_split_subsample.csv")
results_split[, Method := "2. Subsampled split LRT"]

results_hybrid <- fread("sim_data/fig07_hybrid.csv")
results_hybrid[, Method := "3. Subsampled hybrid LRT"]

results <- rbind(results_intersect, results_split, results_hybrid, fill = T) %>% 
  mutate(d = factor(d, levels = c(2, 10, 100, 1000),
                    labels = c("d = 2", "d = 10", "d = 100", "d = 1000"))) %>% 
  filter(n == 1000) %>% 
  as.data.table() 

#########################################################################
##### Plot power at B = 100 with different line types and no points #####
#########################################################################

power_plot_linetype <- results %>%
  filter(theta_star_1 < 0.5 | theta_star_1 > 1) %>% 
  mutate(theta_star_range = factor(theta_star_1 < 0.5, levels = c(T, F),
                                   labels = c("< 0.5", "> 1"))) %>% 
  ggplot(aes(x = theta_star_1, y = reject, col = Method, linetype = Method)) +
  geom_line(data = results %>% 
              filter(theta_star_1 < 0.5, d == "d = 2", theta_star_1 >= 0)) +
  geom_line(data = results %>% 
              filter(theta_star_1 > 1, d == "d = 2", theta_star_1 <= 2)) +
  geom_line(data = results %>% 
              filter(theta_star_1 < 0.5, d == "d = 10", theta_star_1 >= 0)) +
  geom_line(data = results %>% 
              filter(theta_star_1 > 1, d == "d = 10", theta_star_1 <= 2)) +
  geom_line(data = results %>% filter(theta_star_1 < 0.5, d == "d = 100")) +
  geom_line(data = results %>% 
              filter(theta_star_1 > 1, d == "d = 100", theta_star_1 <= 2)) +
  geom_line(data = results %>% filter(theta_star_1 < 0.5, d == "d = 1000")) +
  geom_line(data = results %>% 
              filter(theta_star_1 > 1, d == "d = 1000", theta_star_1 <= 4)) +
  facet_wrap(~ factor(paste(d, theta_star_1 > 1),
                      levels = c("d = 2 FALSE", "d = 2 TRUE",
                                 "d = 10 FALSE", "d = 10 TRUE", 
                                 "d = 100 FALSE", "d = 100 TRUE",
                                 "d = 1000 FALSE", "d = 1000 TRUE"),
                      labels = c("d = 2", "d = 2 ", "d = 10", "d = 10 ",
                                 "d = 100", "d = 100 ", "d = 1000", "d = 1000 ")), 
             scales = "free_x", nrow = 2) +
  labs(title = expression("Power of Tests of"~
                            H[0]*": 0.5 \U2264 ||"*theta*"\U002A|| \U2264 1 vs"~
                            H[1]*": ||"*theta*"\U002A|| < 0.5 or ||"*
                            theta*"\U002A|| > 1"),
       subtitle = "Subsampling at B = 100. n = 1000 observations.",
       y = "Estimated power",
       x = expression(theta[1]*"\U002A, where "~theta*"\U002A"~
                        "="~"("*theta[1]*"\U002A"*", 0, ..., 0)")) +
  scale_color_manual(values = c("black", "blue", "red")) +
  paper_theme +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 35),
        legend.key.width = unit(1.2, "cm"))

###################################################
##### Plot power at B = 100 (journal version) #####
###################################################

power_plot_spaced_journal <- 
  ggplot_gtable(ggplot_build(power_plot_linetype + labs(title = NULL, subtitle = NULL)))

power_plot_spaced_journal$widths[7] <- 0.75*power_plot_spaced_journal$widths[7]
power_plot_spaced_journal$widths[11] <- 1.5*power_plot_spaced_journal$widths[11]
power_plot_spaced_journal$widths[15] <- 0.75*power_plot_spaced_journal$widths[15]
power_plot_spaced_journal$widths[21] <- 2*power_plot_spaced_journal$widths[21]

#####################
##### Save plot #####
#####################

cairo_pdf("sim_plots/figure_07.pdf", height = 5, width = 9)
grid.draw(power_plot_spaced_journal)
dev.off()
