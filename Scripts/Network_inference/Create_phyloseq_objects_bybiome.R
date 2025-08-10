# Import libraries
library(phyloseq)
library(dplyr)        # filter and reformat data frames
library(tibble)       # Needed for converting column to row names
library(microbiome)
library(stringr)
 
 process_biome_data <- function(biome_name, base_folder = "Data"){
   # Create file paths using the biome_experiment_type name
   metadata_file_path <- file.path(base_folder, paste0(biome_name, "/",biome_name, "_metadata.csv"))
   abundance_file_path <- file.path(base_folder, paste0(biome_name, "/",biome_name, "_abundance_table_filtered.csv"))
   taxonomic_file_path <- file.path(base_folder, paste0(biome_name, "/",biome_name, "_taxonomic_table.csv"))
   
   # Load data
   samples_metadata_df <- read.csv(metadata_file_path)
   abundance_df_species <- read.csv(abundance_file_path)
   taxonomic_df_species <- read.csv(taxonomic_file_path)
  
  # Define the row names from the sample column
  samples_metadata_df <- samples_metadata_df %>% 
    tibble::column_to_rownames("Library") 
  
  # Set OTU as rowname
  taxonomic_df_species <- taxonomic_df_species %>% 
    tibble::column_to_rownames("OTU")
  
  abundance_df_species <- abundance_df_species %>% 
    tibble::column_to_rownames("OTU")
  
  # Transform into matrices otu and tax tables
  abund_mat_species <- as.matrix(abundance_df_species)
  tax_mat_species <- as.matrix(taxonomic_df_species)
  
  # Transform to phyloseq objects
  Abund_species = otu_table(abund_mat_species, taxa_are_rows = TRUE)
  Tax_species = tax_table(tax_mat_species)
  samples = sample_data(samples_metadata_df)
  
  phyloseq_file <- phyloseq(Abund_species, Tax_species, samples)

  print(paste0("Phyloseq object created for ", biome_name))
  
  return(phyloseq_file)
}

# List of biome_experiment_types
biome_names <- list.dirs(path = "Data", full.names = FALSE, recursive = FALSE)

# Process each biome
for (biome in biome_names) {
  # Process the data for the biome_experiment_type
  biome_phyloseq <- process_biome_data(biome)
  
  # Construct the file path for saving the phuloseq object
  file_path <- paste0("Data/", biome, "/", biome, "_phyloseqfile_filtered.rds")

  print(paste0("Phyloseq object saved for ", biome))
  
  # Save the biome_experiment_type_phyloseq object
  saveRDS(biome_phyloseq, file = file_path)
}

