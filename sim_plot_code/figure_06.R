# Plot simulated, true, or approximate power of LRTs

library(tidyverse)
library(data.table)
library(gtable)
library(cowplot)
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

# Function to position legend into wasted space.
# Source: https://stackoverflow.com/questions/54438495/shift-legend-into-empty-facets-of-a-faceted-plot-in-ggplot2
shift_legend <- function(p){
  
  # check if p is a valid object
  if(!"gtable" %in% class(p)){
    if("ggplot" %in% class(p)){
      gp <- ggplotGrob(p) # convert to grob
    } else {
      message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
      return(p)
    }
  } else {
    gp <- p
  }
  
  # check for unfilled facet panels
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  if(length(empty.facet.panels) == 0){
    message("There are no unfilled facet panels to shift legend into. Returning original plot.")
    return(p)
  }
  
  # establish extent of unfilled facet panels (including any axis cells in between)
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  empty.facet.panels <- list(min(empty.facet.panels[["t"]]), min(empty.facet.panels[["l"]]),
                             max(empty.facet.panels[["b"]]), max(empty.facet.panels[["r"]]))
  names(empty.facet.panels) <- c("t", "l", "b", "r")
  
  # extract legend & copy over to location of unfilled facet panels
  guide.grob <- which(gp[["layout"]][["name"]] == "guide-box")
  if(length(guide.grob) == 0){
    message("There is no legend present. Returning original plot.")
    return(p)
  }
  gp <- gtable_add_grob(x = gp,
                        grobs = gp[["grobs"]][[guide.grob]],
                        t = empty.facet.panels[["t"]],
                        l = empty.facet.panels[["l"]],
                        b = empty.facet.panels[["b"]],
                        r = empty.facet.panels[["r"]],
                        name = "new-guide-box")
  
  # squash the original guide box's row / column (whichever applicable)
  # & empty its cell
  guide.grob <- gp[["layout"]][guide.grob, ]
  if(guide.grob[["l"]] == guide.grob[["r"]]){
    gp <- gtable_squash_cols(gp, cols = guide.grob[["l"]])
  }
  if(guide.grob[["t"]] == guide.grob[["b"]]){
    gp <- gtable_squash_rows(gp, rows = guide.grob[["t"]])
  }
  gp <- gtable_remove_grobs(gp, "guide-box")
  
  return(gp)
}

# Read in data
usual_lrt <- fread(file = "data/power/02_normal_usual_LRT.csv") %>% 
  mutate(Method = "Classical LRT")

split_lrt <- fread(file = "data/power/03_split_LRT.csv") %>% 
  mutate(Method = "Split LRT",
         power_approx_2 = power)

crossfit_lrt <- fread(file = "data/power/04_crossfit_LRT.csv") %>% 
  mutate(Method = "Cross-fit LRT",
         power_approx_2 = power)

subsample_lrt <- fread(file = "data/power/05_subsample_LRT.csv") %>% 
  mutate(Method = "Subsampling LRT")

# Merge datasets
power_data <- rbind(usual_lrt, split_lrt, crossfit_lrt, subsample_lrt)

# Plot power across LRT methods
power_vs_value <- power_data %>% 
  filter(theta_star <= 0.2) %>% 
  mutate(Method = factor(Method, 
                         levels = c("Classical LRT", "Subsampling LRT",
                                    "Cross-fit LRT", "Split LRT"), 
                         labels = c("Classical LRT", "Subsampling LRT",
                                    "Cross-fit LRT", "Split LRT")),
         d = factor(d, levels = c(1, 2, 10), 
                    labels = c("d = 1", "d = 2", "d = 10"))) %>% 
  ggplot(aes(x = theta_star, y = power, col = Method)) +
  facet_wrap(. ~ d, nrow = 2) +
  geom_line(alpha = 0.8) +
  scale_color_manual(values = c("black", "blue", "red", "orange")) +
  labs(y = "Power",
       x = expression(theta*"*"~"(Repeated in d Dimensions)"),
       title = expression(atop("Power of Tests of H"[0]*":"~theta~"= 0 vs H"[1]*
                                 ":"~theta~"\u2260"~"0",
                               "across Varying True"~theta*"*"))) +
  paper_theme

# Plot power across LRT methods, w/ x axis as distance between true theta and 0
power_vs_norm <- power_data %>% 
  mutate(Method = factor(Method, 
                         levels = c("Classical LRT", "Subsampling LRT",
                                    "Cross-fit LRT", "Split LRT"), 
                         labels = c("Classical LRT", "Subsampling LRT",
                                    "Cross-fit LRT", "Split LRT")),
         d_cat = factor(d, levels = c(1, 2, 10), 
                        labels = c("d = 1", "d = 2", "d = 10")),
         theta_star_sqnorm = theta_star^2 * d) %>% 
  filter(theta_star_sqnorm <= 0.1) %>% 
  ggplot(aes(x = theta_star_sqnorm, y = power, col = Method)) +
  facet_wrap(. ~ d_cat, nrow = 1) +
  geom_line(alpha = 0.8) +
  scale_color_manual(values = c("black", "blue", "red", "orange")) +
  labs(y = "Power",
       x = expression("Squared norm of true mean: ||"*theta*"*"*"||"^{"2"}),
       title = expression(atop("Power of Tests of H"[0]*":"~theta*
                                 "* = 0 vs H"[1]*":"~theta*"* \u2260"~"0",
                               "across Varying True ||"*theta*"*"*"||"^{"2"}))) +
  paper_theme +
  theme(legend.position = "bottom")

power_vs_norm_2 <- power_data %>% 
  mutate(Method = factor(Method, 
                         levels = c("Classical LRT", "Subsampling LRT",
                                    "Cross-fit LRT", "Split LRT"), 
                         labels = c("Classical LRT", "Subsampling LRT",
                                    "Cross-fit LRT", "Split LRT")),
         d_cat = factor(d, levels = c(1, 2, 10), 
                        labels = c("d = 1", "d = 2", "d = 10")),
         theta_star_sqnorm = theta_star^2 * d) %>% 
  filter(theta_star_sqnorm <= 0.1) %>% 
  ggplot(aes(x = theta_star_sqnorm, y = power_approx_2, col = Method)) +
  facet_wrap(. ~ d_cat, nrow = 1) +
  geom_line(alpha = 0.8) +
  scale_color_manual(values = c("black", "blue", "red", "orange")) +
  labs(y = "Power",
       x = expression("Squared norm of true mean: ||"*theta*"*"*"||"^{"2"}),
       title = expression("Power of Tests of H"[0]*":"~theta*
                            "* = 0 vs H"[1]*":"~theta*"* \u2260"~"0"~
                            "across Varying True ||"*theta*"*"*"||"^{"2"}),
       subtitle = expression("Using Normal CDF Approx to Non-central"~
                               chi^2~"CDF")) +
  paper_theme +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 10))

power_vs_norm_2_linetype <- power_data %>% 
  mutate(Method = factor(Method, 
                         levels = c("Classical LRT", "Subsampling LRT",
                                    "Cross-fit LRT", "Split LRT"), 
                         labels = c("Classical LRT", "Subsampling LRT",
                                    "Cross-fit LRT", "Split LRT")),
         d_cat = factor(d, levels = c(1, 2, 10), 
                        labels = c("d = 1", "d = 2", "d = 10")),
         theta_star_sqnorm = theta_star^2 * d) %>% 
  filter(theta_star_sqnorm <= 0.1) %>% 
  ggplot(aes(x = theta_star_sqnorm, y = power_approx_2, col = Method,
             linetype = Method)) +
  facet_wrap(. ~ d_cat, nrow = 1) +
  geom_line(alpha = 0.8) +
  scale_color_manual(values = c("black", "blue", "red", "orange")) +
  scale_linetype_manual(values = c("twodash", "dotdash", "longdash", "solid")) +
  labs(y = "Power",
       x = expression("Squared norm of true mean: ||"*theta*"*"*"||"^{"2"}),
       title = expression("Power of Tests of H"[0]*":"~theta*
                            "* = 0 vs H"[1]*":"~theta*"* \u2260"~"0"~
                            "across Varying True ||"*theta*"*"*"||"^{"2"}),
       subtitle = expression("Using Normal CDF Approx to Non-central"~
                               chi^2~"CDF")) +
  paper_theme +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 10),
        legend.key.width = unit(1.2, "cm"))

######################
##### Save plots #####
######################

ggsave(plot = shift_legend(power_vs_value),
       filename = "plots/power/06_power_vs_value.pdf",
       width = 6, height = 5.4, device = cairo_pdf)

ggsave(plot = power_vs_norm,
       filename = "plots/power/06_power_vs_norm.pdf",
       width = 8, height = 4, device = cairo_pdf)

ggsave(plot = power_vs_norm_2,
       filename = "plots/power/06_power_vs_norm_2.pdf",
       width = 8, height = 4, device = cairo_pdf)

ggsave(plot = power_vs_norm_2 + labs(title = NULL, subtitle = NULL),
       filename = "plots/power/06_power_vs_norm_2_journal.pdf",
       width = 8, height = 3, device = cairo_pdf)

ggsave(plot = power_vs_norm_2_linetype + labs(title = NULL, subtitle = NULL),
       filename = "plots/power/06_power_vs_norm_2_linetype_journal.pdf",
       width = 8, height = 3, device = cairo_pdf)
