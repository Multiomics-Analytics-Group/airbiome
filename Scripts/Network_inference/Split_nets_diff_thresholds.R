library(readr)
library(dplyr)

# Function to load the microNet object
loadNetwork <- function(file_path) {
  net <- readRDS(file_path)
  return(net)
}

# Function to modify and save the microNet object and edge list
modify_and_save_microNet <- function(microNet, thresh, network_type, biome, root_path) {
    # Filter the edge list
    filtered_edges <- microNet$edgelist1 %>% 
      filter(asso >= thresh | asso <= -thresh)
    
    # Skip if no edges
    if (nrow(filtered_edges) == 0) {
      message(paste("No edges for threshold", thresh, "- skipping"))
      return()
    }
    
    # Update the microNet object with the new edge list
    modified_microNet <- microNet
    modified_microNet$edgelist1 <- filtered_edges
    
    # Create folder
    folder_name <- sprintf("%s_%0.2f", network_type, thresh)
    folder_path <- file.path(root_path, biome, folder_name)
    
    if (!dir.exists(folder_path)) {
      dir.create(folder_path)
    }
    
    # Save the RDS file
    rds_file_name <- sprintf("%s_net_%s_%02.0f.rds", biome, network_type, thresh * 100)
    rds_path <- file.path(root_path, biome, folder_name, rds_file_name)
    saveRDS(modified_microNet, file = rds_path)
    
    # Export the network as edge list
    edgelist_file_name <- sprintf("%s_net_%s_%02.0f_edgelist.csv", biome, network_type, thresh * 100)
    edgelist_path <- file.path(root_path, biome, folder_name, edgelist_file_name)
    write.csv(filtered_edges, edgelist_path, row.names = FALSE)
  
    message(paste("microNet object and edge list saved for threshold", thresh))
}

# Base path, network types, and thresholds
root_path <- "Output/Biomes_exptypes"
network_type <- "cclasso"
biomes_extypes <- c("Wastewater_Activated_Sludge_assembly", "Wastewater_Activated_Sludge_metagenomic",
                    "Wastewater_Activated_Sludge_metatranscriptomic", "Wastewater_assembly",
                    "Wastewater_Industrial_wastewater_metagenomic","Wastewater_metagenomic", 
                    "Wastewater_metatranscriptomic", "Wastewater_Water_and_sludge_metagenomic")
  
#c("Wastewater_Water_and_sludge", "Wastewater_Industrial_wastewater", 
#                    "Wastewater_Activated_Sludge", "Wastewater", "assembly",
#                    "metagenomic", "metatranscriptomic")

# biome <- "Wastewater_Water_and_sludge"
thresh_in <- 0.1

# Loop through each biome and threshold
for (biome in biomes_extypes) {
  # Construct file path for the network
  folder_name <- sprintf("%s_%0.1f", network_type, thresh_in)
  net_file_name <- sprintf("%s_net_%s_%02.0f.rds", biome, network_type, thresh_in * 10)
  net_path <- file.path(root_path, biome, folder_name, net_file_name)
  
  # Load 0.1 thres net
  if (file.exists(net_path)) {
    net_0_1 <- loadNetwork(net_path)
    
  } else {
    cat("Network file not found: ", net_path, "\n")
  }
  
  # Apply for each threshold
  thresholds <- seq(0.15, 0.70, by = 0.05)
  
  for (thresh in thresholds){
    # Apply the function to filter out edges and export files
    modify_and_save_microNet(net_0_1, thresh, network_type, biome, root_path)
  }
}
