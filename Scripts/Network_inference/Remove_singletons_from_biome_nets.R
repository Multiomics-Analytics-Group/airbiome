library(igraph)
library(dplyr)

# Define base path, biomes, network types, and thresholds
root_path <- "Output/Biomes_exptypes"
network_type <- "cclasso"
biomes_extypes <- c("Wastewater_Activated_Sludge_assembly", "Wastewater_Activated_Sludge_metagenomic",
                    "Wastewater_Activated_Sludge_metatranscriptomic", "Wastewater_assembly",
                    "Wastewater_Industrial_wastewater_metagenomic","Wastewater_metagenomic", 
                    "Wastewater_metatranscriptomic", "Wastewater_Water_and_sludge_metagenomic")

#  c("Wastewater_Water_and_sludge", "Wastewater_Industrial_wastewater", 
#                    "Wastewater_Activated_Sludge", "Wastewater", "assembly",
#                    "metagenomic", "metatranscriptomic")

thresholds <- seq(0.10, 0.70, by = 0.05)

# Loop through biomes and thresholds
for (biome in biomes_extypes) {
  for (thresh in thresholds) {
    # Construct file path
    folder_name <- sprintf("%s_%0.2f", network_type, thresh)
    folder_path <- file.path(root_path, biome, folder_name)
    net_path <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)
    net_path <- net_path[!grepl("nosinglt", net_path)]
    
    # Check if net_path is empty
    if (length(net_path) == 0) {
      cat("No network files found in folder: ", folder_path, "\n")
      next  # Skip to the next iteration
    }
    
    # Load network if it exists and has more than 3 edges
    if (file.exists(net_path)) {
      net <- read.csv(net_path)
      if (nrow(net) > 3) {
        cat("Network file found and has more than 3 edges: ", net_path, "\n")
        
      } else {
        cat("Network file found but has 3 or fewer edges, skipping: ", net_path, "\n")
        next  # Skip to the next iteration
      }
    } else {
      cat("Network file not found: ", net_path, "\n")
      next  # Skip to the next iteration
    }
    
    # Convert to igraph object
    igraph_net <- graph_from_data_frame(net, directed = FALSE)
    
    # Continuously remove singletons until none are left
    repeat {
      # Calculate degrees
      V(igraph_net)$degree <- degree(igraph_net)
      
      # Identify singletons - nodes with degree < 2
      nodes_to_remove <- V(igraph_net)[degree(igraph_net) < 2]
      
      # Break loop if no singletons
      if (length(nodes_to_remove) == 0) {
        break
      }
      
      # Remove singletons from the graph
      igraph_net <- delete_vertices(igraph_net, nodes_to_remove)
    }
    
    # Conver net to edge list
    net_nosinglet_edgelist <- igraph::as_data_frame(igraph_net, what = "edges")
    
    # Change column names from 'from' and 'to' to 'v1' and 'v2'
    colnames(net_nosinglet_edgelist) <- c("v1", "v2", "asso", "diss", "adja")
    
    # Export network without singletons as and edgelist
    output_filename <- sprintf("%s_nosinglt_%s_%02.0f_edgelist.csv",biome, network_type, thresh * 100)
    output_filepath <- file.path(folder_path, output_filename)
    write.csv(net_nosinglet_edgelist, output_filepath, row.names = FALSE)
    
  }
}
