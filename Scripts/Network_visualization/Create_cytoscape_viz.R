library(RCy3)
library(igraph)
library(dplyr)
library(optparse) 
library(paletteer) 

# Set seed for reproducibility
set.seed(123)

# Define a function to visualize the network in Cytoscape
visualize_network_cytoscape <- function(net_path_arg, weight_col_arg, collection_name_arg) {
  # --- Input Validation and Loading ---
  if (!file.exists(net_path_arg)) {
    stop(paste0("Error: Network file not found at ", net_path_arg))
  }

  net <- read.csv(net_path_arg)

  if (nrow(net) <= 3) {
    stop(paste0("Error: Network has 3 or fewer edges, skipping visualization: ", net_path_arg))
  }

  # Ensure the weight column exists in the loaded data
  if (!weight_col_arg %in% colnames(net)) {
    stop(paste0("Error: Weight column '", weight_col_arg, "' not found in the network file."))
  }

  cat(paste0("Loading and processing network: ", net_path_arg, "\n"))

  # Add 'sign' attribute to the data frame if 'asso' is the weight column, BEFORE converting to igraph.
  # This ensures Cytoscape receives the 'sign' attribute upon network creation.
  if (weight_col_arg == "asso") {
    if (!"asso" %in% colnames(net)) {
      stop("Error: 'asso' column not found in data frame for signed coloring.")
    }
    net$sign <- sign(net$asso)
    cat("  'sign' column added to data frame for signed coloring in Cytoscape.\n")
  }

  # Convert to igraph object
  # Ensure the graph is created with the specified weight_col_arg as an edge attribute
  igraph_net <- graph_from_data_frame(net, directed = FALSE)

  # --- Node and Edge Attribute Calculations ---
  # Add node degree as an attribute
  V(igraph_net)$degree <- degree(igraph_net)
  min_degree <- min(V(igraph_net)$degree)
  max_degree <- max(V(igraph_net)$degree)

  # Calculate clusters (communities)
  clusters <- tryCatch({
    cluster_louvain(igraph_net) 
  }, error = function(e) {
    message("Warning: Louvain community detection failed, likely due to disconnected graph. Assigning default clusters.")
    # Assign each node to its own cluster if community detection fails
    list(membership = seq_len(vcount(igraph_net)))
  })
  V(igraph_net)$cluster <- clusters$membership

  # Generate a color palette for clusters
  num_clusters <- length(unique(V(igraph_net)$cluster))
  base_palette <- paletteer_d("ggsci::category20_d3", 20)
  base_palette <- substr(base_palette, 1, 7) # Ensure hex codes are clean

  # Extend the color palette by repeating it if necessary
  repeat_times <- ceiling(num_clusters / length(base_palette))
  color_palette <- rep(base_palette, times = repeat_times)[1:num_clusters]

  # Create a mapping between cluster number and color
  cluster_color_mapping <- setNames(color_palette, unique(V(igraph_net)$cluster))

  # --- Cytoscape Visualization ---
  # IMPORTANT: This script assumes Cytoscape Desktop application is running and CyREST API is enabled.

  # Derive network title from the file path
  network_title <- basename(net_path_arg)
  network_title <- tools::file_path_sans_ext(network_title) # Remove .csv extension

  cat(paste0("Creating network '", network_title, "' in Cytoscape collection '", collection_name_arg, "'\n"))
  RCy3::createNetworkFromIgraph(igraph_net, title = network_title, collection = collection_name_arg)

  # Get network SUID for applying styles
  networkSuid <- RCy3::getNetworkSuid()

  # Create a new visual style or reuse if exists
  style_name <- paste0("Style_", network_title, "_", weight_col_arg)
  if (!(style_name %in% getVisualStyleNames())) {
    createVisualStyle(style_name)
  }
  setVisualStyle(style_name, network = networkSuid)
  
  cat(paste0("Applying style '", style_name, "' to network.\n"))

  # Set default node shape
  RCy3::setNodeShapeDefault('ellipse', style.name = style_name)

  # Set node size mapping based on degree
  RCy3::setNodeSizeMapping('degree', c(min_degree, round(max_degree * 0.2), max_degree), c(10, 20, 60), 
                            network = networkSuid, style.name = style_name)

  # Set node color mapping based on cluster
  RCy3::setNodeColorMapping('cluster', as.numeric(names(cluster_color_mapping)), mapping.type = "d",
                            cluster_color_mapping, network = networkSuid, style.name = style_name)

  # Set node label mapping to the 'name' attribute (which comes from v1/v2 columns in input CSV)
  RCy3::setNodeLabelMapping('name', network = networkSuid, style.name = style_name)
  cat("Node labels set to node 'name' attribute.\n")

  # --- Conditional Edge Color Mapping based on Weight Column ---
  # Define colors for positive and negative edges (relevant for 'asso' type)
  positive_color <- "#ac1214"  # Red for positive
  negative_color <- "#0d1fad"  # Blue for negative
  default_edge_color <- "#808080" # Grey for non-signed or default edges

  if (weight_col_arg == "asso") {
    RCy3::setEdgeColorMapping('sign', c(-1, 1), c(negative_color, positive_color),
                              network = networkSuid, style.name = style_name)
    cat("  Signed edge color mapping applied.\n")
    
    # Optionally, also set edge width based on absolute 'asso' strength
    edge_weights_abs <- abs(E(igraph_net)$asso)
    min_weight_abs <- min(edge_weights_abs)
    max_weight_abs <- max(edge_weights_abs)
    RCy3::setEdgeLineWidthMapping(weight_col_arg, 
                                 c(min_weight_abs, max_weight_abs), 
                                 c(1, 8), # Min and max line width
                                 network = networkSuid, style.name = style_name, mapping.type = "c")
  } else {
    # For 'adja' or 'diss', do not consider sign; apply a default color or a gradient for positive values
    cat(paste0("Applying default edge color for '", weight_col_arg, "' (no sign considered).\n"))
    # If the edge attribute exists and has numeric values, map to a continuous color scale
    if (weight_col_arg %in% names(E(igraph_net)$attr) && is.numeric(E(igraph_net)$attr[[weight_col_arg]])) {
      edge_values <- E(igraph_net)$attr[[weight_col_arg]]
      min_val <- min(edge_values)
      max_val <- max(edge_values)
      
      # Use a single-hue gradient for positive weights (e.g., from light blue to dark blue)
      RCy3::setEdgeColorMapping(weight_col_arg, 
                                c(min_val, max_val), 
                                c("#e0f2f7", "#0d47a1"), # Light blue to dark blue gradient
                                network = networkSuid, style.name = style_name, mapping.type = "c")
      
      RCy3::setEdgeLineWidthMapping(weight_col_arg, 
                                 c(min_val, max_val), 
                                 c(1, 8), # Min and max line width
                                 network = networkSuid, style.name = style_name, mapping.type = "c")
    } else {
      # Fallback to a single default color if the attribute is not numeric or not found
      RCy3::setEdgeColorDefault(default_edge_color, style.name = style_name)
      RCy3::setEdgeLineWidthDefault(2, style.name = style_name) # Default width
    }
  }

  cat("Network visualization complete in Cytoscape.\n")
}

# --- Command-line Argument Parsing ---
option_list <- list(
  make_option(c("-n", "--net_path"), type="character", 
              help="Full path to the input network CSV file.",
              metavar="file_path"),
  make_option(c("-w", "--weight_column"), type="character", default="adja",
              help="Column to use as edge weight ('adja', 'asso', or 'diss'). Affects coloring and size mapping. [default %default]",
              metavar="weight_col"),
  make_option(c("-c", "--collection_name"), type="character", default="MicrobialNetworks",
              help="Name of the Cytoscape collection to add the network to. [default %default]",
              metavar="collection_name")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# --- Main Execution ---
if (is.null(opt$net_path)) {
  stop("Error: --net_path argument is required. Please specify the path to your network CSV file.", call. = FALSE)
}

# Run the visualization function with parsed arguments
tryCatch({
  visualize_network_cytoscape(opt$net_path, opt$weight_column, opt$collection_name)
}, error = function(e) {
  message(paste0("An error occurred during network visualization: ", e$message))
})