library(igraph)
library(dplyr)
library(optparse)

# ----- Main -----
# Define command-line options
option_list <- list(
  make_option(c("-r", "--root_path"), type="character", default="Output/Network_inference/NOPGE_sp_01", 
              help="Root directory where filtered network files are located and where cleaned networks will be saved.", 
              metavar="directory"),
  make_option(c("-t", "--network_type"), type="character", default="cclasso",
              help="Type of network (e.g., 'cclasso' or 'spring'). Used for folder naming.", 
              metavar="type"),
  make_option(c("--start_thresh"), type="numeric", default=0.10, 
              help="Starting threshold for the loop to find network files (e.g., 0.10).", 
              metavar="number"),
  make_option(c("--end_thresh"), type="numeric", default=0.95, 
              help="Ending threshold for the loop to find network files (e.g., 0.95).", 
              metavar="number"),
  make_option(c("--step_thresh"), type="numeric", default=0.05, 
              help="Step size for the threshold loop (e.g., 0.05).", 
              metavar="number")
)

# Parse the command-line options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Use the parsed options
root_path <- opt$root_path
network_type <- opt$network_type
start_thresh <- opt$start_thresh
end_thresh <- opt$end_thresh
step_thresh <- opt$step_thresh

# Define thresholds sequence based on parsed arguments
thresholds <- seq(start_thresh, end_thresh, by = step_thresh)

# Loop through thresholds to find and process network files
for (thresh in thresholds) {
  # Construct the expected folder path for network files at the current threshold
  folder_name <- sprintf("%s_%0.2f", network_type, thresh)
  folder_path <- file.path(root_path, folder_name)
  
  # List CSV files in the folder, excluding those already marked as 'nosinglt' (no singletons)
  net_path_candidates <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)
  net_path_candidates <- net_path_candidates[!grepl("nosinglt", net_path_candidates)]
  
  # Check if any network file was found for this threshold and folder
  if (length(net_path_candidates) == 0) {
    cat(paste0("No uncleaned network CSV files found in folder: ", folder_path, "\n"))
    next # Skip to the next iteration if no files are found
  }
  
  # Assuming there's only one relevant CSV file per folder (or picking the first one)
  net_path <- net_path_candidates[1] 

  # Load network if it exists and has more than 3 edges (a common check to ensure meaningful networks)
  if (file.exists(net_path)) {
    net <- read.csv(net_path)
    if (nrow(net) > 3) {
      cat(paste0("Processing network file: ", net_path, " (edges: ", nrow(net), ")\n"))
    } else {
      cat(paste0("Network file found but has 3 or fewer edges, skipping: ", net_path, "\n"))
      next # Skip to the next iteration if not enough edges
    }
  } else {
    cat(paste0("Network file not found at expected path: ", net_path, "\n"))
    next # Skip to the next iteration if the file doesn't exist
  }

  # Convert the data frame to an igraph object for graph operations
  igraph_net <- graph_from_data_frame(net, directed = FALSE)

  # Continuously remove singletons (nodes with degree < 2) until none are left
  # A "singleton" here is a node that is not connected to at least one other node.
  repeat {
    # Calculate the degree (number of connections) for each vertex (node) in the graph
    V(igraph_net)$degree <- degree(igraph_net)
    
    # Identify nodes that have a degree less than 2 (i.e., isolated nodes or those with only one connection)
    nodes_to_remove <- V(igraph_net)[degree(igraph_net) < 2]
    
    # If no more singletons are found, exit the loop
    if (length(nodes_to_remove) == 0) {
      break
    }
    
    # Remove the identified singleton nodes from the graph
    igraph_net <- delete_vertices(igraph_net, nodes_to_remove)
  }

  # Convert the cleaned igraph object back to an edge list data frame
  net_nosinglet_edgelist <- igraph::as_data_frame(igraph_net, what = "edges")

  # Merge the cleaned edgelist (from igraph) with the preserved columns
  final_edgelist <- net_nosinglet_edgelist %>%
    semi_join(net %>% select(v1, v2, asso, diss, adja), by = c("from" = "v1", "to" = "v2")) %>%
    # Ensure correct column names if they were not already 'from', 'to' in igraph output
    select(v1 = from, v2 = to, asso, diss, adja) # Explicitly select and rename

  # Export the network without singletons as a new edgelist CSV file
  output_filename <- sprintf("%s_nosinglt_%02.0f_edgelist.csv", network_type, thresh * 100)
  output_filepath <- file.path(folder_path, output_filename)
  write.csv(final_edgelist, output_filepath, row.names = FALSE)
  
  cat(paste0("Cleaned network saved to: ", output_filepath, "\n"))
}
