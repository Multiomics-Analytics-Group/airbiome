# Import libraries
library(phyloseq)
library(dplyr)       # filter and reformat data frames
library(tibble)      # Needed for converting column to row names
library(microbiome)  # Potentially useful for further analysis, though not directly used in this snippet
library(stringr)     # Potentially useful for string manipulation, though not directly used in this snippet

# Define a function to create a phyloseq object from a given folder
create_phyloseq_from_folder <- function(folder_path) {
  
  folder_name <- basename(folder_path)
  print(paste0("Processing data for folder: ", folder_name))
  
  # --- Identify file paths based on patterns ---
  all_files_in_folder <- list.files(path = folder_path, full.names = TRUE)
  
  # Find metadata file
  metadata_files <- all_files_in_folder[grepl("metadata", all_files_in_folder, ignore.case = TRUE)]
  if (length(metadata_files) == 1) {
    metadata_file_path <- metadata_files[1]
  } else {
    stop(paste0("Error: Expected exactly one file containing 'metadata' in its name in ", folder_path, ", but found ", length(metadata_files), "."))
  }
  
  # Find abundance file
  abundance_files <- all_files_in_folder[grepl("abundance_table_filtered", all_files_in_folder, ignore.case = TRUE)]
  if (length(metadata_files) == 1) {
    abundance_file_path <- abundance_files[1]
  } else {
    stop(paste0("Error: Expected exactly one file containing 'abundance_table_filtered' in its name in ", folder_path, ", but found ", length(abundance_files), "."))
  }
  
  # Find taxonomic file
  taxonomic_files <- all_files_in_folder[grepl("taxonomic_table_filtered", all_files_in_folder, ignore.case = TRUE)]
  if (length(taxonomic_files) == 1) {
    taxonomic_file_path <- taxonomic_files[1]
  } else {
    stop(paste0("Error: Expected exactly one file containing 'taxonomic_table_filtered' in its name in ", folder_path, ", but found ", length(taxonomic_files), "."))
  }
  
  print(paste0("Identified metadata file: ", basename(metadata_file_path)))
  print(paste0("Identified abundance file: ", basename(abundance_file_path)))
  print(paste0("Identified taxonomic file: ", basename(taxonomic_file_path)))
  
  # --- Load data ---
  samples_metadata_df <- read.csv(metadata_file_path)
  abundance_df <- read.csv(abundance_file_path)
  taxonomic_df <- read.csv(taxonomic_file_path)
  
  # --- Define the row names ---
  # For samples_metadata_df, "Library" is the column to be set as row names
  if ("Library" %in% colnames(samples_metadata_df)) {
    samples_metadata_df <- samples_metadata_df %>%
      tibble::column_to_rownames("Library")
  } else {
    stop("Error: 'Library' column not found in metadata.csv. Cannot set row names.")
  }
  
  # Determine the ID column name dynamically by looking for a column ending with "_ID"
  id_cols <- colnames(taxonomic_df)[grepl("_ID$", colnames(taxonomic_df))]
  
  if (length(id_cols) == 0) { # Check if no ID column is found at all
    stop("Error: No column ending with '_ID' found in the taxonomic table.")
  }
  id_column_name <- id_cols[1] 
  
  # For abundance_df and taxonomic_df, use the dynamically determined ID column
  if (id_column_name %in% colnames(taxonomic_df)) {
    taxonomic_df <- taxonomic_df %>%
      tibble::column_to_rownames(id_column_name)
  } else {
    stop(paste0("Error: '", id_column_name, "' column not found in the taxonomic table. Cannot set row names."))
  }
  
  if (id_column_name %in% colnames(abundance_df)) {
    abundance_df <- abundance_df %>%
      tibble::column_to_rownames(id_column_name)
  } else {
    stop(paste0("Error: '", id_column_name, "' column not found in the abundance table. Cannot set row names."))
  }
  
  # --- Transform into matrices for phyloseq ---
  abund_mat <- as.matrix(abundance_df)
  tax_mat <- as.matrix(taxonomic_df)
  
  # --- Transform to phyloseq objects ---
  Abund = otu_table(abund_mat, taxa_are_rows = TRUE)
  Tax = tax_table(tax_mat)
  samples = sample_data(samples_metadata_df)
  
  # Create the phyloseq object
  phyloseq_file <- phyloseq(Abund, Tax, samples)
  
  print(paste0("Phyloseq object created for folder: ", folder_name))
  
  # --- Save the phyloseq object ---
  # Construct the file path for saving the phyloseq object
  file_path <- file.path(folder_path, paste0(folder_name, "_phyloseqfile_filtered", ".rds"))
  
  # Save the phyloseq object
  saveRDS(phyloseq_file, file = file_path)
  
  print(paste0("Phyloseq object saved to: ", file_path))
  
  return(phyloseq_file) # Optionally return the phyloseq object
}

# --- Example Calls ---
# Call the function with the species folder_path
pge_sp_folder_path <- "Data/PGE_new/Processed_data/PGE_Species" 
phyloseq_object_pge_sp <- create_phyloseq_from_folder(pge_sp_folder_path)

nopge_sp_folder_path <- "Data/NOPGE_new/Processed_data/NOPGE_Species" 
phyloseq_object_nopge_sp <- create_phyloseq_from_folder(nopge_sp_folder_path)

# Call the function with the genus folder_path
pge_genus_folder_path <- "Data/PGE_new/Processed_data/PGE_Genus" 
phyloseq_object_pge_gen <- create_phyloseq_from_folder(pge_genus_folder_path)

nopge_genus_folder_path <- "Data/NOPGE_new/Processed_data/NOPGE_Genus" 
phyloseq_object_nopge_gen <- create_phyloseq_from_folder(nopge_genus_folder_path)

