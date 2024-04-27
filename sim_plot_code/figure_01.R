# Plot showing regions of LRT coverage

library(tidyverse)
library(data.table)
# library(devtools)
# devtools::install_github("cmartin/ggConvexHull")
library(ggConvexHull)

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
grid_data_usual <- 
  fread(file = "sim_data/fig01_classical_LRT.csv") %>% 
  mutate(sim = factor(sim, levels = 1:6,
                      labels = c("Sim 1", "Sim 2", "Sim 3",
                                 "Sim 4", "Sim 5", "Sim 6")),
         Method = "Classical LRT") %>%
  dplyr::select(coord1, coord2, Method, covered, sim) %>% 
  as.data.table()

grid_data_split <- 
  fread(file = "sim_data/fig01_split_LRT.csv") %>% 
  mutate(sim = factor(sim, levels = 1:6,
                      labels = c("Sim 1", "Sim 2", "Sim 3",
                                 "Sim 4", "Sim 5", "Sim 6"))) %>%
  mutate(Method = "Split LRT",
         covered = split_LRT_cov) %>% 
  dplyr::select(coord1, coord2, Method, covered, sim) %>% 
  as.data.table()

grid_data_crossfit <- 
  fread(file = "sim_data/fig01_crossfit_LRT.csv") %>% 
  mutate(sim = factor(sim, levels = 1:6,
                      labels = c("Sim 1", "Sim 2", "Sim 3",
                                 "Sim 4", "Sim 5", "Sim 6"))) %>%
  mutate(Method = "Cross-fit LRT",
         covered = crossfit_cov) %>% 
  dplyr::select(coord1, coord2, Method, covered, sim) %>% 
  as.data.table()

grid_data_subsample <- 
  fread(file = "sim_data/fig01_subsample_LRT.csv") %>% 
  mutate(sim = factor(sim, levels = 1:6,
                      labels = c("Sim 1", "Sim 2", "Sim 3",
                                 "Sim 4", "Sim 5", "Sim 6")),
         Method = "Subsampling LRT") %>%
  dplyr::select(coord1, coord2, Method, covered, sim) %>% 
  as.data.table()

grid_data <- rbind(grid_data_usual, grid_data_split,
                   grid_data_crossfit, grid_data_subsample) %>% 
  mutate(Method = factor(Method, 
                         levels = c("Classical LRT", "Subsampling LRT", 
                                    "Cross-fit LRT", "Split LRT"),
                         labels = c("Classical LRT", "Subsampling LRT", 
                                    "Cross-fit LRT", "Split LRT")))

########################
##### Create plots #####
########################

grid_plot <- grid_data %>%
  ggplot(aes(x = coord1, y = coord2, col = Method)) +
  geom_convexhull(alpha = 0, aes(fill = Method)) +
  facet_wrap(. ~ sim) +
  labs(x = "Coordinate 1",
       y = "Coordinate 2",
       title = "Coverage Regions of Multivariate Normal LRT Methods at d = 2",
       subtitle = expression("Fixed"~paste(theta, "*")~"= (0, 0),"~
                               "Fixed Data (n = 1000), Random Subsampling")) +
  scale_color_manual(values = c("black", "blue", "red", "orange")) +
  paper_theme +
  theme(legend.text = element_text(margin = margin(t = 2, b = 10, unit = "pt")),
        legend.key.size = unit(1, "cm"),
        aspect.ratio = 1,
        axis.text = element_text(size = 9))

######################
##### Save plots #####
######################

ggsave(plot = grid_plot + labs(title = NULL, subtitle = NULL),
       filename = "sim_plots/figure_01.pdf",
       width = 8.5, height = 4.5)
