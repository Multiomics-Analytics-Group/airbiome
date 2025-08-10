library(readr)
library(dplyr)
library(optparse)

# Function to load the microNet object
loadNetwork <- function(file_path) {
  net <- readRDS(file_path)
  return(net)
}

# Function to filter and save rds network object and CSV per threshold
modify_and_save_microNet <- function(microNet, thresh, network_type, root_path, filter_col) {
  # Dynamic filtering based on the specified column
  if (filter_col == "asso") {
    # For 'asso' column, filter for values >= thresh OR <= -thresh
    filtered_edges <- microNet$edgelist %>%
      filter(!!sym(filter_col) >= thresh | !!sym(filter_col) <= -thresh)
  } else {
    # For 'adja' or 'diss' (or any other column), filter for values >= thresh
    filtered_edges <- microNet$edgelist %>%
      filter(!!sym(filter_col) >= thresh)
  }

  # Skip if no edges
  if (nrow(filtered_edges) == 0) {
    message(paste("No edges for threshold", thresh, "- skipping"))
    return()
  }

  # Update the microNet object with filtered edges
  filtered_net <- microNet
  filtered_net$edgelist <- filtered_edges

  # Create output folder
  folder_path <- file.path(root_path, sprintf("%s_%0.2f", network_type, thresh))
  dir.create(folder_path, recursive = TRUE, showWarnings = FALSE)

  # File names
  base_name <- sprintf("%s_net_%02.0f", network_type, thresh * 100)
  rds_path <- file.path(folder_path, paste0(base_name, ".rds"))
  csv_path <- file.path(folder_path, paste0(base_name, "_edgelist.csv"))

  # Save files
  saveRDS(filtered_net, rds_path)
  write.csv(filtered_edges, csv_path, row.names = FALSE)

  message(sprintf("microNet object and edge list (csv) saved for threshold %.2f", thresh))
}

# ----- Main -----
# Define command-line options
option_list <- list(
  make_option(c("-r", "--root_path"), type="character", default="Output/Network_inference/NOPGE_sp_01", 
              help="Root directory for output", metavar="directory"),
  make_option(c("-n", "--net_file_name"), type="character", default="NOPGE_sp_net_cclasso_01.rds", 
              help="Network filename", metavar="file"),
  make_option(c("-t", "--net_type"), type="character", default="cclasso", 
              help="Type of network", metavar="string"),
  make_option(c("-i", "--thresh_input"), type="numeric", default=0.1, 
              help="Input threshold", metavar="numeric"),
  make_option(c("-c", "--filter_col"), type="character", default="adja", 
              help="Column to filter on: 'adja', 'asso', or 'diss'", metavar="string")
)

# Parse the command-line options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Use the parsed options
root_path <- opt$root_path
net_file_name <- opt$net_file_name
network_type <- opt$net_type
thresh_input <- opt$thresh_input
filter_col <- opt$filter_col

# Construct the full network path
folder_name <- sprintf("%s_%0.1f", network_type, thresh_input)
net_path <- file.path(root_path, folder_name, net_file_name)

# Load network
if (!file.exists(net_path)) {
  stop("Network file not found: ", net_path, call. = FALSE)
}
net_01 <- loadNetwork(net_path)

# Apply for each threshold
thresholds <- seq(0.15, 0.90, by = 0.05)
for (thresh in thresholds) {
  modify_and_save_microNet(net_01, thresh, network_type, root_path, filter_col)
}
