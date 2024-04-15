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

results_intersect <- fread("data/RIPR/norm_interval/01_intersect.csv")
results_intersect[, Method := "1. Intersection"]

results_split <- fread("data/RIPR/norm_interval/02_split_subsample.csv")
results_split[, Method := "2. Subsampled split LRT"]

results_hybrid <- fread("data/RIPR/norm_interval/03_hybrid.csv")
results_hybrid[, Method := "3. Subsampled hybrid LRT"]

results <- rbind(results_intersect, results_split, results_hybrid, fill = T) %>% 
  mutate(d = factor(d, levels = c(2, 10, 100, 1000),
                    labels = c("d = 2", "d = 10", "d = 100", "d = 1000"))) %>% 
  filter(n == 1000) %>% 
  as.data.table() 

#################################
##### Plot power at B = 100 #####
#################################

power_plot <- results %>%
  filter(theta_star_1 < 0.5 | theta_star_1 > 1) %>% 
  mutate(theta_star_range = factor(theta_star_1 < 0.5, levels = c(T, F),
                                   labels = c("< 0.5", "> 1"))) %>% 
  ggplot(aes(x = theta_star_1, y = reject, col = Method)) +
  geom_point(alpha = 0.5,
             data = results %>% 
               filter(theta_star_1 < 0.5, d == "d = 2", theta_star_1 >= 0)) +
  geom_line(data = results %>% 
              filter(theta_star_1 < 0.5, d == "d = 2", theta_star_1 >= 0)) +
  geom_point(alpha = 0.5,
             data = results %>% 
               filter(theta_star_1 > 1, d == "d = 2", theta_star_1 <= 2)) +
  geom_line(data = results %>% 
              filter(theta_star_1 > 1, d == "d = 2", theta_star_1 <= 2)) +
  geom_point(alpha = 0.5,
             data = results %>% 
               filter(theta_star_1 < 0.5, d == "d = 10", theta_star_1 >= 0)) +
  geom_line(data = results %>% 
              filter(theta_star_1 < 0.5, d == "d = 10", theta_star_1 >= 0)) +
  geom_point(alpha = 0.5,
             data = results %>% 
               filter(theta_star_1 > 1, d == "d = 10", theta_star_1 <= 2)) +
  geom_line(data = results %>% 
              filter(theta_star_1 > 1, d == "d = 10", theta_star_1 <= 2)) +
  geom_point(alpha = 0.5,
             data = results %>% filter(theta_star_1 < 0.5, d == "d = 100")) +
  geom_line(data = results %>% filter(theta_star_1 < 0.5, d == "d = 100")) +
  geom_point(alpha = 0.5,
             data = results %>% 
               filter(theta_star_1 > 1, d == "d = 100", theta_star_1 <= 2)) +
  geom_line(data = results %>% 
              filter(theta_star_1 > 1, d == "d = 100", theta_star_1 <= 2)) +
  geom_point(alpha = 0.5,
             data = results %>% filter(theta_star_1 < 0.5, d == "d = 1000")) +
  geom_line(data = results %>% filter(theta_star_1 < 0.5, d == "d = 1000")) +
  geom_point(alpha = 0.5,
             data = results %>% 
               filter(theta_star_1 > 1, d == "d = 1000", theta_star_1 <= 4)) +
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
        axis.text.x = element_text(size = 10, angle = 35))

power_plot_spaced <- ggplot_gtable(ggplot_build(power_plot))

power_plot_spaced$widths[7] <- 0.75*power_plot_spaced$widths[7]
power_plot_spaced$widths[11] <- 1.5*power_plot_spaced$widths[11]
power_plot_spaced$widths[15] <- 0.75*power_plot_spaced$widths[15]
power_plot_spaced$widths[21] <- 2*power_plot_spaced$widths[21]

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

power_plot_spaced <- ggplot_gtable(ggplot_build(power_plot_linetype))

power_plot_spaced$widths[7] <- 0.75*power_plot_spaced$widths[7]
power_plot_spaced$widths[11] <- 1.5*power_plot_spaced$widths[11]
power_plot_spaced$widths[15] <- 0.75*power_plot_spaced$widths[15]
power_plot_spaced$widths[21] <- 2*power_plot_spaced$widths[21]

###################################################
##### Plot power at B = 100 (journal version) #####
###################################################

power_plot_spaced_journal <- 
  ggplot_gtable(ggplot_build(power_plot_linetype + labs(title = NULL, subtitle = NULL)))

power_plot_spaced_journal$widths[7] <- 0.75*power_plot_spaced_journal$widths[7]
power_plot_spaced_journal$widths[11] <- 1.5*power_plot_spaced_journal$widths[11]
power_plot_spaced_journal$widths[15] <- 0.75*power_plot_spaced_journal$widths[15]
power_plot_spaced_journal$widths[21] <- 2*power_plot_spaced_journal$widths[21]

#############################
##### Plot type I error #####
#############################

validity_plot <- results %>% 
  filter(theta_star_1 >= 0.5, theta_star_1 <= 1) %>% 
  ggplot(aes(x = theta_star_1, y = reject, col = Method)) +
  geom_point(alpha = 0.5) +
  geom_line() +
  geom_hline(yintercept = 0.10, lty = "dashed") +
  scale_y_continuous(limits = c(0, 0.2)) +
  facet_wrap(~ d) +
  labs(title = expression(atop("False Positive Rate of Tests of"~
                                 H[0]*": 0.5 \U2264 ||"*theta*"\U002A|| \U2264 1 vs",
                               H[1]*": ||"*theta*"\U002A|| < 0.5 or ||"*
                                 theta*"\U002A|| > 1")),
       subtitle = "Subsampling at B = 100. n = 1000 observations.",
       y = "False Positive Rate",
       x = expression(theta[1]*"\U002A, where "~theta*"\U002A"~
                        "="~"("*theta[1]*"\U002A"*", 0, ..., 0)")) +
  scale_color_manual(values = c("black", "blue", "red")) +
  paper_theme +
  theme(legend.position = "bottom")

#############################################################
##### Plot case 1/2/3 proportions for hybrid split/RIPR #####
#############################################################

case_prop_plot <- results %>%
  filter(Method == "3. Subsampled hybrid LRT",
         !(d == "d = 1000" &
             as.character(theta_star_1) %in% c(seq(1.91, 1.99, by = 0.01),
                                               seq(2.01, 2.05, by = 0.01)))) %>% 
  pivot_longer(cols = c("mean_case_1_prop", "mean_case_2_prop",
                        "mean_case_3_prop"),
               names_to = "Case", values_to = "mean_prop") %>% 
  mutate(Case = factor(Case, 
                       levels = c("mean_case_1_prop", "mean_case_2_prop",
                                  "mean_case_3_prop"),
                       labels = c(1, 2, 3))) %>% 
  ggplot(aes(x = theta_star_1, y = mean_prop, col = Case)) +
  geom_point(alpha = 0.8) +
  geom_line(alpha = 0.8) +
  facet_wrap(~ d, scales = "free") +
  scale_color_manual(
    values = c("darkgreen", "darkorange", "purple"),
    labels = c(expression("(1)  ||"*bar(Y)["1,b"]*"|| < 0.5   "), 
               expression("(2)  0.5 \u2264 ||"*bar(Y)["1,b"]*"|| \u2264 1   "),
               expression("(3)  ||"*bar(Y)["1,b"]*"|| > 1"))) +
  labs(title = "Proportions of three cases for hybrid LRT",
       subtitle = "Subsampling at B = 100. n = 1000 observations.",
       y = "Proportion of simulations",
       x = expression(theta[1]*"\U002A, where "~theta*"\U002A"~
                        "="~"("*theta[1]*"\U002A"*", 0, ..., 0)")) +
  paper_theme +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 35))

# Journal version (without titles)

case_prop_plot_journal <- case_prop_plot + 
  labs(title = NULL,
       subtitle = NULL)

######################
##### Save plots #####
######################

cairo_pdf("plots/RIPR/norm_interval/04_power.pdf", height = 5.5, width = 9)
grid.draw(power_plot_spaced)
dev.off()

cairo_pdf("plots/RIPR/norm_interval/04_power_journal.pdf", height = 5, width = 9)
grid.draw(power_plot_spaced_journal)
dev.off()

ggsave(plot = validity_plot, 
       filename = "plots/RIPR/norm_interval/04_validity.pdf",
       width = 9, height = 5, device = cairo_pdf)

ggsave(plot = case_prop_plot, 
       filename = "plots/RIPR/norm_interval/04_case_prop.pdf",
       width = 8, height = 5.5, device = cairo_pdf)

ggsave(plot = case_prop_plot_journal, 
       filename = "plots/RIPR/norm_interval/04_case_prop_journal.pdf",
       width = 8, height = 5, device = cairo_pdf)
