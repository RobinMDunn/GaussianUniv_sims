# GaussianUniv Repository Overview

This repository contains code to replicate the results of [Gaussian universal likelihood ratio testing](https://arxiv.org/abs/2104.14676) by [Robin Dunn](https://robinmdunn.github.io/), [Aaditya Ramdas](http://www.stat.cmu.edu/~aramdas/), [Sivaraman Balakrishnan](https://www.stat.cmu.edu/~siva/), and [Larry Wasserman](https://www.stat.cmu.edu/~larry/). 

## Folder contents

- [sim_code](sim_code): Code for the paper's simulations. Scripts are labeled by the figure for which they simulate data. Each R script saves the simulation output to [sim_data](sim_data). The [combine_datasets.R](sim_code/combine_datasets.R) script
provides a template to combine the output from multiple simulations into a
single dataset.
- [sim_data](sim_data): Output of simulations from [sim_code](sim_code).
- [sim_plot_code](sim_plot_code): Code to reproduce the paper's plots and tables. Reads in data from [sim_data](sim_data) and outputs plots to [sim_plots](sim_plots).
- [sim_plots](sim_plots): Plots from the paper. The plots are the output of the scripts in [sim_plot_code](sim_plot_code).
