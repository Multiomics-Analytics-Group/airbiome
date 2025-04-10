# Use user-level library path
user_lib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE)
.libPaths(user_lib)

# Import libraries
# List of required packages
required_packages <- c("devtools", "BiocManager", "optparse")

# Install missing packages from CRAN
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos = 'https://cloud.r-project.org')

# Install missing packages from Bioconductor
BiocManager::install(c("phyloseq", "microbiome"))

# Install NetCoMi dependencies
devtools::install_github("zdk123/SpiecEasi")
devtools::install_github("GraceYoon/SPRING")

# Install NetCoMi
devtools::install_github("stefpeschel/NetCoMi", 
                         repos = c("https://cloud.r-project.org/",
                                   BiocManager::repositories()))

# Add NetCoMi to the list of required packages
required_packages <- c(required_packages, "NetCoMi", "phyloseq", "microbiome")

# Load the packages
lapply(required_packages, require, character.only = TRUE)

# Set seed to make results reproducible
set.seed(123)

# Define function to create and save networks
constructAndSaveNetwork <- function(networkType, data1, data2, thresh, taxRank="Species", seed=123, 
                                    base_folder, file_name_prefix, biome_name) {
    # Common parameters
    common_params <- list(
      taxRank = taxRank,
      normMethod = "none",
      dissFunc = "signed",
      verbose = 3,
      seed = seed
    )
    
    # Network-specific parameters
    if (networkType == "spring") {
      network_params <- c(list(
        measure = "spring",
        measurePar = list(nlambda = 20, rep.num = 20, thresh = thresh),
        zeroMethod = "none",
        sparsMethod = "none"
      ), common_params)
    } else if (networkType == "cclasso") {
      network_params <- c(list(
        measure = "cclasso",
        zeroMethod = "pseudoZO",
        sparsMethod = "threshold",
        thresh = thresh
      ), common_params)
    } else {
      stop("Unsupported network type")
    }
    
    # Construct the network
    net <- do.call(netConstruct, c(list(data = data1, data2 = data2), network_params))
    
    # Check if the network has edges
    if (is.null(net$edgelist1) || nrow(net$edgelist1) == 0) {
      message("Network ", networkType, " at threshold ", thresh, " for biome ", 
              biome_name, " has no edges. Skipping.")
      return(NULL)
    }
    
    # Save the RDS file
    rds_file_name <- file.path(base_folder, paste0(biome_name, file_name_prefix, ".rds"))
    saveRDS(net, file = rds_file_name)
    
    # Export the networks as edge lists
    net_df <- as.data.frame(do.call(cbind, net$edgelist1))
    csv_file_name <- file.path(base_folder, paste0(biome_name, file_name_prefix, "_edgelist.csv"))
    write.csv(net_df, csv_file_name, row.names = FALSE)
}

# Define the options
option_list <- list(
  make_option(c("-t", "--threshold"), type="numeric", default=10, 
              help="Threshold value", metavar="number"),
  make_option(c("-n", "--network"), type="character", default="cclasso",
              help="Network type ('cclasso' or 'spring')", metavar="type"),
  make_option(c("-b", "--biome"), type="character", default="PGE_and_NOPGE",
              help="Biome", metavar="biome")
)

# Parse the options
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate network type and biome
valid_network_types <- c("cclasso", "spring")
valid_biomes <- c("PGE_and_NOPGE")

if (!opt$network %in% valid_network_types) {
  stop("Network type must be 'cclasso' or 'spring'")
}

if (!opt$biome %in% valid_biomes) {
  stop("Biome must be one of the specified options")
}

# Load both group datasets
phyloseq_group1 <- readRDS("Data/PGE_phyloseqfile_filtered.rds")
phyloseq_group2 <- readRDS("Data/NO_PGE_phyloseqfile_filtered.rds")

# Create file_path
file_output_path <- paste0("Output/", opt$biome)

# Check if the directory exists, if not create it
if (!dir.exists(file_output_path)) {
  dir.create(file_output_path, recursive = TRUE)
}

# Create network
net_path <- file.path(file_output_path, sprintf("%s_%0.1f", opt$network, opt$threshold))
network_prefix <- sprintf("_net_%s_%02.0f", opt$network, opt$threshold * 10)

# Check if the directory exists, if not create it
if (!dir.exists(net_path)) {
  dir.create(net_path, recursive = TRUE)
}

constructAndSaveNetwork(opt$network, phyloseq_group1, phyloseq_group2, opt$threshold, base_folder=net_path, 
                        file_name_prefix=network_prefix, biome_name=opt$biome)