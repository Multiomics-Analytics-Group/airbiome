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
                      centrLCC = TRUE,
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
                      #removeSingle = TRUE,
                      sPathNorm = TRUE)
  return(props)
}

# Function to summarize network properties
summarizeProperties <- function(net, props, file_path) {
  summary_text <- capture.output({
    
    # First network
    edgelist <- net$edgelist1
    degree_data <- as.data.frame(props$centralities$degree1)
    colnames(degree_data) <- "degree"
    
    num_nodes <- length(unique(c(edgelist$v1, edgelist$v2)))
    num_edges <- nrow(edgelist)
    avg_degree <- degree_data %>%
      dplyr::filter(degree != 0) %>%
      dplyr::summarize(avg_degree = mean(degree, na.rm = TRUE)) %>%
      dplyr::pull(avg_degree)
    
    cat("Number of Nodes: ", num_nodes, "\n")
    cat("Number of Edges: ", num_edges, "\n")
    cat("Average Degree (excluding singletons): ", avg_degree, "\n\n")
    
    # Combined summary
    cat("=== Network Properties Summary ===\n")
    print(summary(props, clusterLCC = TRUE, numbNodes = 10L, digits = 3L))
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
root_path <- "/Users/asaru/Documents/DTU/MoNA/Airbiome/Output/Best_nets"
network_type <- "cclasso"
biome <- "PGE_and_NOPGE_comb_gen"
thresh <- 0.7

# Loop through each network type and threshold
folder_name <- "PGE_and_NOPGE_combined_gen"
net_file_name <- "PGE_NOPGE_gen_net_cclasso_07.rds"
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

} else {
  cat("Network file not found: ", net_path, "\n")
}

quant_comp <- netCompare(props, permTest = FALSE, 
                          verbose = FALSE, seed = 222)

saveRDS(quant_comp, file = "/Users/asaru/Documents/DTU/MoNA/Airbiome/Output/quant_comp_airbiome_cclasso_05.rds")

saveRDS(props, file = "/Users/asaru/Documents/DTU/MoNA/Airbiome/Output/desc_comp_airbiome_cclasso_05.rds")


diff_net_PGE <- diffnet(net, diffMethod = "fisherTest", adjust = "lfdr")

