library(RCy3)
library(igraph)
library(dplyr)
library(optparse)
library(paletteer)

# Set seed for reproducibility
set.seed(123)

# Define a function to visualize a single network in Cytoscape
visualize_single_network <- function(net_path, weight_column, biome_name, threshold, network_type) {
  cat(paste0("\n--- Processing Network: Biome='", biome_name, "', Threshold='", threshold, "', Type='", network_type, "' ---\n"))
  
  # --- Input Validation and Loading ---
  if (!file.exists(net_path)) {
    cat(paste0("Error: Network file not found at ", net_path, ". Skipping.\n"))
    return(NULL)
  }

  net_df <- read.csv(net_path)

  if (nrow(net_df) <= 3) {
    cat(paste0("Error: Network has 3 or fewer edges (", nrow(net_df), " edges) at ", net_path, ". Skipping visualization.\n"))
    return(NULL)
  }

  if (!weight_column %in% colnames(net_df)) {
    stop(paste0("Error: Weight column '", weight_column, "' not found in the network file: ", net_path))
  }

  cat(paste0("  Loading network: ", net_path, "\n"))

  # Add 'sign' attribute to the data frame if 'asso' is the weight column, BEFORE converting to igraph.
  # This ensures Cytoscape receives the 'sign' attribute upon network creation.
  if (weight_column == "asso") {
    if (!"asso" %in% colnames(net_df)) {
      stop("Error: 'asso' column not found in data frame for signed coloring.")
    }
    net_df$sign <- sign(net_df$asso)
    cat("  'sign' column added to data frame for signed coloring.\n")
  }

  # Convert to igraph object
  igraph_net <- graph_from_data_frame(net_df, directed = FALSE)
  cat("  igraph object created. Nodes: ", vcount(igraph_net), ", Edges: ", ecount(igraph_net), "\n")

  # --- Node and Edge Attribute Calculations ---
  # Add node degree as an attribute
  V(igraph_net)$degree <- degree(igraph_net)
  min_degree <- min(V(igraph_net)$degree)
  max_degree <- max(V(igraph_net)$degree)
  cat(paste0("  Node degrees calculated. Min degree: ", min_degree, ", Max degree: ", max_degree, ".\n"))

  # Calculate clusters (communities)
  clusters <- tryCatch({
    cluster_louvain(igraph_net) 
  }, error = function(e) {
    message("Warning: Louvain community detection failed, likely due to disconnected graph or no edges. Assigning each node to its own cluster.")
    list(membership = seq_len(vcount(igraph_net))) # Assign each node to its own cluster if community detection fails
  })
  V(igraph_net)$cluster <- clusters$membership
  cat("  Clusters calculated and assigned to nodes. Number of clusters: ", length(unique(V(igraph_net)$cluster)), "\n")

  # Generate a color palette for clusters
  num_clusters <- length(unique(V(igraph_net)$cluster))
  base_palette <- paletteer_d("ggsci::category20_d3", 20)
  base_palette <- substr(base_palette, 1, 7) # Ensure hex codes are clean

  repeat_times <- ceiling(num_clusters / length(base_palette))
  color_palette <- rep(base_palette, times = repeat_times)[1:num_clusters]

  cluster_color_mapping <- setNames(color_palette, unique(V(igraph_net)$cluster))
  cat("  Cluster color mapping created.\n")

  # --- Cytoscape Visualization ---
  # IMPORTANT: This script assumes Cytoscape Desktop application is running and CyREST API is enabled.
  # Define collection name using biome for better organization
  collection_name <- paste0(collection_name_arg, "_", biome_name)

  # Derive network title for Cytoscape
  network_title <- paste0(network_type, "_", sprintf("%02.0f", threshold * 100), "_", biome_name)
  cat(paste0("  Creating network '", network_title, "' in Cytoscape collection '", collection_name, "'\n"))
  RCy3::createNetworkFromIgraph(igraph_net, title = network_title, collection = collection_name)

  # Get network SUID for applying styles
  networkSuid <- RCy3::getNetworkSuid()
  cat(paste0("  Network SUID: ", networkSuid, "\n"))

  # Create a new visual style or reuse if exists
  style_name <- paste0("Style_", network_type, "_", sprintf("%02.0f", threshold * 100), "_", biome_name, "_", weight_column)
  if (!(style_name %in% getVisualStyleNames())) {
    createVisualStyle(style_name)
    cat(paste0("  Created new visual style: '", style_name, "'\n"))
  } else {
    cat(paste0("  Reusing existing visual style: '", style_name, "'\n"))
  }
  setVisualStyle(style_name, network = networkSuid)
  cat(paste0("  Applied style '", style_name, "' to network.\n"))

  # Set default node shape
  RCy3::setNodeShapeDefault('ellipse', style.name = style_name)
  cat("  Default node shape set to 'ellipse'.\n")

  # Set node size mapping based on degree
  cat("  Setting node size mapping based on degree...\n")
  if (min_degree == max_degree) {
    RCy3::setNodeSizeDefault(20, style.name = style_name)
    cat("  Warning: All nodes have the same degree. Node size set to a default of 20.\n")
  } else {
    tryCatch({
      RCy3::setNodeSizeMapping('degree', c(min_degree, round(max_degree * 0.2), max_degree), c(10, 20, 60),
                                network = networkSuid, style.name = style_name)
      cat("  Node size mapping (degree) applied successfully.\n")
    }, error = function(e) {
      message(paste0("  ERROR applying node size mapping (degree): ", e$message))
      message("  Setting default node size due to mapping error.")
      RCy3::setNodeSizeDefault(20, style.name = style_name)
    })
  }

  # Set node color mapping based on cluster
  cat("  Setting node color mapping based on cluster...\n")
  RCy3::setNodeColorMapping('cluster', as.numeric(names(cluster_color_mapping)), mapping.type = "d",
                            cluster_color_mapping, network = networkSuid, style.name = style_name)
  cat("  Node color mapping (cluster) applied.\n")

  # Set node label mapping to the 'name' attribute
  cat("  Setting node label mapping to node 'name' attribute...\n")
  RCy3::setNodeLabelMapping('name', network = networkSuid, style.name = style_name)
  cat("  Node labels set.\n")

  # --- Conditional Edge Color Mapping based on Weight Column ---
  positive_color <- "#ac1214"  # Red for positive
  negative_color <- "#0d1fad"  # Blue for negative
  default_edge_color <- "#808080" # Grey for non-signed or default edges

  if (weight_column == "asso") {
    cat("  'asso' weight column selected. Applying signed edge coloring.\n")
    RCy3::setEdgeColorMapping('sign', c(-1, 1), c(negative_color, positive_color),
                              network = networkSuid, style.name = style_name)
    cat("  Signed edge color mapping applied.\n")

    edge_weights_abs <- abs(igraph::edge_attr(igraph_net, "asso"))
    min_weight_abs <- min(edge_weights_abs)
    max_weight_abs <- max(edge_weights_abs)

    cat("  Applying edge width mapping based on absolute 'asso' strength.\n")
    RCy3::setEdgeLineWidthMapping("asso",
                                 c(min_weight_abs, max_weight_abs),
                                 c(1, 8), # Min and max line width
                                 network = networkSuid, style.name = style_name, mapping.type = "c")
    cat("  Edge width mapping applied.\n")
  } else {
    cat(paste0("  '", weight_column, "' weight column selected. Applying default edge coloring.\n"))
    if (weight_column %in% igraph::edge_attr_names(igraph_net) && is.numeric(igraph::edge_attr(igraph_net, weight_column))) {
      edge_values <- igraph::edge_attr(igraph_net, weight_column)
      min_val <- min(edge_values)
      max_val <- max(edge_values)

      cat("  Applying continuous edge color mapping based on numeric weight values.\n")
      RCy3::setEdgeColorMapping(weight_column,
                                c(min_val, max_val),
                                c("#e0f2f7", "#0d47a1"), # Light blue to dark blue gradient
                                network = networkSuid, style.name = style_name, mapping.type = "c")

      cat("  Applying continuous edge width mapping based on numeric weight values.\n")
      RCy3::setEdgeLineWidthMapping(weight_column,
                                 c(min_val, max_val),
                                 c(1, 8), # Min and max line width
                                 network = networkSuid, style.name = style_name, mapping.type = "c")
    } else {
      cat("  Weight attribute not numeric or found, falling back to default edge style.\n")
      RCy3::setEdgeColorDefault(default_edge_color, style.name = style_name)
      RCy3::setEdgeLineWidthDefault(2, style.name = style_name) # Default width
    }
  }

  cat("--- Network visualization complete in Cytoscape ---\n")
}

# --- Command-line Argument Parsing ---
option_list <- list(
  make_option(c("-r", "--root_path"), type="character",
              help="Base directory containing biome subfolders (e.g., '/path/to/Output/Network_inference/Indiv_biomes').",
              metavar="root_path"),
  make_option(c("-t", "--network_type"), type="character", default="cclasso",
              help="Type of network to process (e.g., 'cclasso', 'spring'). [default %default]",
              metavar="network_type"),
  make_option(c("-b", "--biomes"), type="character",
              help="Comma-separated list of biome names to process (e.g., 'PGE,NOPGE').",
              metavar="biomes"),
  make_option(c("-s", "--thresholds"), type="character",
              help="Comma-separated list of thresholds (e.g., '0.10,0.90,0.30') or a range 'start:end:step' (e.g., '0.1:0.9:0.05').",
              metavar="thresholds"),
  make_option(c("-w", "--weight_column"), type="character", default="adja",
              help="Column to use as edge weight ('adja', 'asso', or 'diss'). Affects coloring and size mapping. [default %default]",
              metavar="weight_col"),
  make_option(c("-c", "--collection_name_prefix"), type="character", default="MicrobialNetworks",
              help="Prefix for Cytoscape collection names (e.g., 'Airbiome_Nets'). Each biome will get its own collection.",
              metavar="collection_name_prefix")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# --- Validate Required Arguments ---
if (is.null(opt$root_path)) {
  stop("Error: --root_path argument is required.", call. = FALSE)
}
if (is.null(opt$biomes)) {
  stop("Error: --biomes argument is required (comma-separated list).", call. = FALSE)
}
if (is.null(opt$thresholds)) {
  stop("Error: --thresholds argument is required (comma-separated list or range).", call. = FALSE)
}

# --- Parse Biomes and Thresholds ---
biomes_list <- unlist(strsplit(opt$biomes, ","))

thresholds_num <- NULL
if (grepl(":", opt$thresholds)) { # Check for range format 'start:end:step'
  parts <- as.numeric(unlist(strsplit(opt$thresholds, ":")))
  if (length(parts) == 3) {
    thresholds_num <- seq(parts[1], parts[2], by = parts[3])
  } else {
    stop("Error: Invalid threshold range format. Use 'start:end:step' or 'val1,val2,val3'.", call. = FALSE)
  }
} else { # Assume comma-separated list
  thresholds_num <- as.numeric(unlist(strsplit(opt$thresholds, ",")))
  if (any(is.na(thresholds_num))) {
    stop("Error: Invalid threshold list. Ensure all values are numeric and comma-separated.", call. = FALSE)
  }
}

# --- Main Loop to Process and Visualize Networks ---
cat("\n--- Starting Batch Network Visualization ---\n")
cat(paste0("Root Path: ", opt$root_path, "\n"))
cat(paste0("Network Type: ", opt$network_type, "\n"))
cat(paste0("Biomes: ", paste(biomes_list, collapse=", "), "\n"))
cat(paste0("Thresholds: ", paste(thresholds_num, collapse=", "), "\n"))
cat(paste0("Weight Column: ", opt$weight_column, "\n"))
cat(paste0("Collection Prefix: ", opt$collection_name_prefix, "\n"))


for (biome in biomes_list) {
  biome_path <- file.path(opt$root_path, biome)
  if (!dir.exists(biome_path)) {
    cat(paste0("Warning: Biome directory not found: ", biome_path, ". Skipping this biome.\n"))
    next
  }

  for (thresh in thresholds_num) {
    # Construct folder name (e.g., cclasso_0.50)
    folder_name <- sprintf("%s_%0.2f", opt$network_type, thresh)
    folder_path <- file.path(biome_path, folder_name)
    
    # Locate the network file (assuming it's the nosinglt CSV within the folder)
    net_files_in_folder <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)
    net_path_to_use <- net_files_in_folder[grepl("nosinglt", net_files_in_folder)]

    if (length(net_path_to_use) == 0) {
      cat(paste0("Warning: No 'nosinglt' network CSV found in ", folder_path, ". Skipping this network.\n"))
      next
    } else if (length(net_path_to_use) > 1) {
      cat(paste0("Warning: Multiple 'nosinglt' network CSVs found in ", folder_path, ". Using the first one: ", basename(net_path_to_use[1]), ".\n"))
      net_path_to_use <- net_path_to_use[1]
    } else {
      net_path_to_use <- net_path_to_use[1]
    }
    
    # Call the single network visualization function
    tryCatch({
      visualize_single_network(net_path_to_use, opt$weight_column, biome, thresh, opt$network_type)
    }, error = function(e) {
      message(paste0("An error occurred while visualizing network for Biome '", biome, "', Threshold '", thresh, "': ", e$message))
    })
  }
}

cat("\n--- All requested networks processed. ---\n")
