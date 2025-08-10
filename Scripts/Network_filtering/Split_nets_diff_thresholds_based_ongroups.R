library(readr)
library(dplyr)

# Load the microNet object
loadNetwork <- function(file_path) {
  net <- readRDS(file_path)
  return(net)
}

modify_and_save_microNet <- function(microNet, thresh, network_type, biome, root_path) {
  # ----- Combined -----
  combined_microNet <- microNet
  combined_microNet$edgelist1 <- microNet$edgelist1 %>%
    filter(adja >= thresh)
  combined_microNet$edgelist2 <- microNet$edgelist2 %>%
    filter(adja >= thresh)

  combined_folder <- file.path(root_path, biome, sprintf("%s_combined_%0.2f", network_type, thresh))
  dir.create(combined_folder, recursive = TRUE, showWarnings = FALSE)
  combined_rds <- file.path(combined_folder, sprintf("%s_net_%s_combined_%02.0f.rds", biome, network_type, thresh * 100))
  saveRDS(combined_microNet, combined_rds)
  message(paste("Saved network rds object for PGE and NOPGE at threshold", thresh))

  # ----- PGE -----
  filtered_edges1 <- combined_microNet$edgelist1
  if (nrow(filtered_edges1) > 0) {
    pge_net <- list(edgelist1 = filtered_edges1)
    pge_folder <- file.path(root_path, "PGE", sprintf("%s_%0.2f", network_type, thresh))
    dir.create(pge_folder, recursive = TRUE, showWarnings = FALSE)

    pge_rds <- file.path(pge_folder, sprintf("PGE_net_%s_%02.0f.rds", network_type, thresh * 100))
    pge_csv <- file.path(pge_folder, sprintf("PGE_net_%s_%02.0f_edgelist.csv", network_type, thresh * 100))
    saveRDS(pge_net, pge_rds)
    write.csv(filtered_edges1, pge_csv, row.names = FALSE)
    message(paste("Saved network and edgelist for PGE at threshold", thresh))

  } else {
    message(paste("No PGE edges for threshold", thresh))
  }

  # ----- NOPGE -----
  filtered_edges2 <- combined_microNet$edgelist2
  if (nrow(filtered_edges2) > 0) {
    nopge_net <- list(edgelist2 = filtered_edges2)
    nopge_folder <- file.path(root_path, "NOPGE", sprintf("%s_%0.2f", network_type, thresh))
    dir.create(nopge_folder, recursive = TRUE, showWarnings = FALSE)

    nopge_rds <- file.path(nopge_folder, sprintf("NOPGE_net_%s_%02.0f.rds", network_type, thresh * 100))
    nopge_csv <- file.path(nopge_folder, sprintf("NOPGE_net_%s_%02.0f_edgelist.csv", network_type, thresh * 100))
    saveRDS(nopge_net, nopge_rds)
    write.csv(filtered_edges2, nopge_csv, row.names = FALSE)
    message(paste("Saved network and edgelist for NOPGE at threshold", thresh))

  } else {
    message(paste("No NOPGE edges for threshold", thresh))
  }

  message(paste("Finished threshold", thresh))
}

# Base path, network types, and thresholds
root_path <- "Output/Network_inference/PGE_and_NOPGE_05"
network_type <- "cclasso"
biome <- "PGE_and_NOPGE"
thresh_in <- 0.5

# Construct file path for the network
folder_name <- sprintf("%s_%0.1f", network_type, thresh_in)
net_file_name <- sprintf("%s_net_%s_%02.0f.rds", biome, network_type, thresh_in * 10)
net_path <- file.path(root_path, folder_name, net_file_name)

# Load 0.5 thres net
if (file.exists(net_path)) {
  net_05 <- loadNetwork(net_path)
  
} else {
  cat("Network file not found: ", net_path, "\n")
}

# Apply for each threshold
thresholds <- seq(0.50, 0.90, by = 0.05)

for (thresh in thresholds){
  # Apply the function to filter out edges and export files
  modify_and_save_microNet(net_05, thresh, network_type, biome, root_path)
}
