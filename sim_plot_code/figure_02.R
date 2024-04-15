# Plot simulated and approximate analytical expressions for 
# subsampling as B goes to infinity.
# Panels considering several choices of d and n.

suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(latex2exp))

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
results <- fread("data/radius/08_ddim_limit_results.csv")

# Plot of simulated expectation and analytical expression values
expectation_plot <- results %>% 
  filter(d %in% c(1, 20), n %in% c(10, 50, 200)) %>% 
  mutate(d_and_n = factor(paste0("d = ", d, ", n = ", n),
                          levels = c("d = 1, n = 10", "d = 1, n = 50",   
                                     "d = 1, n = 200", "d = 20, n = 10",
                                     "d = 20, n = 50", "d = 20, n = 200"),
                          labels = c("d = 1, n = 10", "d = 1, n = 50",   
                                     "d = 1, n = 200", "d = 20, n = 10",
                                     "d = 20, n = 50", "d = 20, n = 200"))) %>% 
  ggplot(aes(x = c_val)) +
  facet_wrap(~ d_and_n, scales = "free") +
  geom_line(data = . %>% filter(method == "Analytical\nexpression"),
            aes(y = value, col = "Analytical\nexpression\n"),
            alpha = 0.7, na.rm = T) +
  geom_point(data = . %>% filter(method == "Simulated\nexpectation"),
             aes(y = value, col = "Simulated\nexpectation\n"),
             alpha = 0.7) +
  geom_hline(yintercept = 10, lty = "dashed", col = "black") +
  labs(title = expression(atop("Simulated Expectation vs Analytical Expression", 
                               "for Limiting Subsampling Test Statistic as"~
                                 B~symbol('\256')~infinity)),
       x = TeX("$c$, where $\\theta = c\\textbf{1}$"),
       y = "Value",
       col = "Method") +
  scale_color_manual(values = c("Simulated\nexpectation\n" = "black",
                                "Analytical\nexpression\n" = "red")) +
  paper_theme +
  theme(axis.text = element_text(size = 10))

# Plot of simulated expectation and analytical expression values,
# with single line and single point for legend
expectation_plot_new_legend <- results %>% 
  filter(d %in% c(1, 20), n %in% c(10, 50, 200)) %>% 
  mutate(d_and_n = factor(paste0("d = ", d, ", n = ", n),
                          levels = c("d = 1, n = 10", "d = 1, n = 50",   
                                     "d = 1, n = 200", "d = 20, n = 10",
                                     "d = 20, n = 50", "d = 20, n = 200"),
                          labels = c("d = 1, n = 10", "d = 1, n = 50",   
                                     "d = 1, n = 200", "d = 20, n = 10",
                                     "d = 20, n = 50", "d = 20, n = 200"))) %>% 
  ggplot(aes(x = c_val)) +
  facet_wrap(~ d_and_n, scales = "free") +
  geom_line(data = . %>% filter(method == "Analytical\nexpression"),
            aes(y = value, linetype = "Analytical\nexpression\n"),
            alpha = 0.7, na.rm = T, col = "red") +
  geom_point(data = . %>% filter(method == "Simulated\nexpectation"),
             aes(y = value, shape = "Simulated\nexpectation\n"),
             alpha = 0.7, col = "black") +
  geom_hline(yintercept = 10, lty = "dashed", col = "black") +
  labs(title = expression(atop("Simulated Expectation vs Analytical Expression", 
                               "for Limiting Subsampling Test Statistic as"~
                                 B~symbol('\256')~infinity)),
       x = TeX("$c$, where $\\theta = c\\textbf{1}$"),
       y = "Value",
       col = "Method") +
  paper_theme +
  scale_linetype_discrete(name = NULL) +
  scale_shape_discrete(name = NULL) +
  theme(axis.text = element_text(size = 10))

######################
##### Save plots #####
######################

ggsave(plot = expectation_plot,
       filename = "plots/radius/08_expect_analytical_panel.pdf",
       width = 9, height = 5)

ggsave(plot = expectation_plot + labs(title = NULL, subtitle = NULL),
       filename = "plots/radius/08_expect_analytical_panel_journal.pdf",
       width = 9, height = 4)

ggsave(plot = expectation_plot_new_legend + labs(title = NULL, subtitle = NULL),
       filename = "plots/radius/08_expect_analytical_new_legend_journal.pdf",
       width = 9, height = 4)
