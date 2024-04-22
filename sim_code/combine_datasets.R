# Merge all .csv files with a given file name prefix into a single .csv file.
# Change start_name (line 8) and file name (line 52) for different data.

# Read in libraries
library(dplyr)
library(data.table)

# Set common file name start for all files to merge
start_name <- "^fig01_classical_LRT_sim_"

# Get all file names with given start_name
filenames <- list.files(path = "sim_data", pattern = start_name, full.names = T)

# Read in first file.
file_1 <- data.table::fread(file = filenames[1])

# Check if dataset includes mean_test_stat column.
# If yes, read in datasets individually, setting this column to numeric.
# (Some large values may have been stored as characters.)
# If no, read in all data with map_df.
if("mean_test_stat" %in% colnames(file_1)) {

  file_list <- vector("list", length(filenames))

  for(i in 1:length(filenames)) {
    file_list[[i]] <- data.table::fread(file = filenames[i])
    file_list[[i]]$mean_test_stat <- as.numeric(file_list[[i]]$mean_test_stat)
  }

  all_data <- do.call("rbind", file_list)

} else {

  all_data <- list.files(path = "sim_data",
                         pattern = start_name, full.names = T) %>%
    purrr::map_df(~data.table::fread(.))

}

# Save new data frame
data.table::fwrite(all_data, file = "sim_data/fig01_classical_LRT.csv")
