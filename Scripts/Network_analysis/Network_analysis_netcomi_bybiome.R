# Import libraries
library(NetCoMi)
library(dplyr)

# Set seed to make results reproducible
set.seed(123)

# Function to load networks
loadNetwork <- function(file_path) {
  net <- readRDS(file_path)
  return(net)
}

# Function to analyze networks 
analyzeNetwork <- function(net) {
  props <- netAnalyze(net, 
                      # Centrality
                      centrLCC = FALSE,
                      weightDeg = TRUE,
                      normDeg = FALSE,
                      
                      # Hubs
                      hubPar = "degree",
                      hubQuant = 0.9,
                      lnormFit = TRUE,
                      
                      # Clustering
                      clustMethod = "cluster_fast_greedy",
                      weightClustCoef = TRUE,
                      
                      # Others
                      avDissIgnoreInf = TRUE,
                      removeSingle = TRUE,
                      sPathNorm = TRUE)
  return(props)
}

# Function to summarize network properties for two edgelists and a joint props object
summarizeProperties <- function(net, props, file_path) {
  summary_text <- capture.output({
    
    # First network
    edgelist1 <- net$edgelist1
    degree_data1 <- as.data.frame(props$centralities$degree1)
    colnames(degree_data1) <- "degree"
    
    num_nodes1 <- length(unique(c(edgelist1$v1, edgelist1$v2)))
    num_edges1 <- nrow(edgelist1)
    avg_degree1 <- degree_data1 %>%
      dplyr::filter(degree != 0) %>%
      dplyr::summarize(avg_degree = mean(degree, na.rm = TRUE)) %>%
      dplyr::pull(avg_degree)
    
    cat("=== Network 1 ===\n")
    cat("Number of Nodes: ", num_nodes1, "\n")
    cat("Number of Edges: ", num_edges1, "\n")
    cat("Average Degree (excluding singletons): ", avg_degree1, "\n\n")
    
    # Second network
    edgelist2 <- net$edgelist2
    degree_data2 <- as.data.frame(props$centralities$degree2)
    colnames(degree_data2) <- "degree"
    
    num_nodes2 <- length(unique(c(edgelist2$v1, edgelist2$v2)))
    num_edges2 <- nrow(edgelist2)
    avg_degree2 <- degree_data2 %>%
      dplyr::filter(degree != 0) %>%
      dplyr::summarize(avg_degree = mean(degree, na.rm = TRUE)) %>%
      dplyr::pull(avg_degree)
    
    cat("=== Network 2 ===\n")
    cat("Number of Nodes: ", num_nodes2, "\n")
    cat("Number of Edges: ", num_edges2, "\n")
    cat("Average Degree (excluding singletons): ", avg_degree2, "\n\n")
    
    # Combined summary
    cat("=== Combined Network Properties Summary ===\n")
    print(summary(props, numbNodes = 10L, digits = 3L))
  })
  
  # Save to file
  writeLines(summary_text, con = file_path)
}

# Function to plot networks
plotNetwork <- function(props, file_path) {
  png(file_path, width=900, height=700)
  
  plot(props, layout = "spring", sameLayout = TRUE, layoutGroup = "union", 
       rmSingles = "inboth", repulsion = 0.8, nodeColor = "cluster", 
       edgeTranspLow = 0, edgeTranspHigh = 40, nodeSize = "degree", posCol = "#ac1214", 
       negCol = "#0d1fad", labels = FALSE, cexNodes = 2, nodeSizeSpread = 2.6, 
       hubBorderCol = "snow4")
  
  legend(-1, 1, title = "Estimated associations:", legend = c("+","-"), 
         col = c("#ac1214","#0d1fad"), inset = 0.02, cex = 1.4, lty = 1, lwd = 3.5, 
         bty = "n", horiz = TRUE, y.intersp = 0.7)
  
  dev.off()
}

# Base path, network types, and thresholds
root_path <- "/Users/asaru/Documents/DTU/MoNA/Airbiome/Output/Network_inference/PGE_and_NOPGE_05"
network_type <- "cclasso"
biome <- "PGE_and_NOPGE"
thresh <- 0.5

# Loop through each network type and threshold
folder_name <- sprintf("%s_%0.1f", network_type, thresh)
net_file_name <- sprintf("%s_net_%s_%02.0f.rds", biome, network_type, thresh * 10)
net_path <- file.path(root_path, folder_name, net_file_name)

# Load, analyze, and summarize properties of the network if the file exists
if (file.exists(net_path)) {
  net <- loadNetwork(net_path)
  props <- analyzeNetwork(net)
  
  # Summarize properties
  summary_file_name <- sprintf("summ_props_%s_%s_%0.1f.txt", biome, 
                               network_type, thresh)
  summary_file_path <- file.path(root_path, folder_name, 
                                 summary_file_name)
  summarizeProperties(net, props, summary_file_path)
  
  # Plot network
  plot_file_name <- sprintf("%s_%s_net_gen_deg_%0.1f.png", biome, network_type, thresh)
  plot_file_path <- file.path(root_path, folder_name, plot_file_name)
  plotNetwork(props, plot_file_path)
} else {
  cat("Network file not found: ", net_path, "\n")
}

quant_comp <- netCompare(props, permTest = FALSE, 
                          verbose = FALSE, seed = 222)

saveRDS(quant_comp, file = "/Users/asaru/Documents/DTU/MoNA/Airbiome/Output/quant_comp_airbiome_cclasso_05.rds")

saveRDS(props, file = "/Users/asaru/Documents/DTU/MoNA/Airbiome/Output/desc_comp_airbiome_cclasso_05.rds")


diff_net_PGE <- diffnet(net, diffMethod = "fisherTest", adjust = "lfdr")

