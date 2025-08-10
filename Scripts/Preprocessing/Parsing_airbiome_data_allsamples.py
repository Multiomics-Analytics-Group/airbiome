import pandas as pd
from collections import defaultdict
from pathlib import Path
import argparse

# Import logging utilities from utils.py
from utils import (
    get_logger,
    load_config,
)


def filt_otus_by_relfreq_and_prev(
    abundance_df, min_rel_freq=0.00005, min_prevalence=0.25
):
    """
    Filters OTUs (Operational Taxonomic Units) from an abundance DataFrame
    based on relative frequency and prevalence thresholds.

    Args:
        abundance_df (pd.DataFrame): DataFrame where index is OTU and columns are samples.
                                     Values are read counts.
        min_rel_freq (float): Minimum relative frequency an OTU must have in a sample
                              to be considered 'present'.
        min_prevalence (float): Minimum proportion of samples an OTU must be 'present' in
                                to be retained.

    Returns:
        pd.DataFrame: A new DataFrame containing only the OTUs that meet the filtering criteria.
                      Returns an empty DataFrame if the input is empty or lacks numeric data.
    """
    if abundance_df.empty:
        logger.info("No OTUs to filter. Returning empty DataFrame.")
        return pd.DataFrame()

    # Ensure all columns are numeric for sum operation
    numeric_cols = abundance_df.select_dtypes(include=["number"]).columns
    if numeric_cols.empty:
        logger.info(
            "No numeric columns found for abundance calculation. Returning empty DataFrame."
        )
        return pd.DataFrame()

    # Calculate relative abundance only for numeric columns
    # Divides each column by its sum to get relative frequencies for each sample
    rel_abundance = abundance_df[numeric_cols].div(
        abundance_df[numeric_cols].sum(axis=0), axis=1
    )

    # Create a boolean matrix indicating where relative abundance meets the minimum frequency
    presence_matrix = rel_abundance >= min_rel_freq
    # Calculate prevalence: sum of 'True' (present) divided by total number of samples
    prevalence = presence_matrix.sum(axis=1) / presence_matrix.shape[1]
    # Get the index (OTU IDs) of OTUs that meet the prevalence threshold
    keep_otus = prevalence[prevalence >= min_prevalence].index
    # Return the original abundance DataFrame, filtered to keep only the selected OTUs
    return abundance_df.loc[keep_otus]


def get_genus_and_species(name):
    """
    Extracts genus and species from a full taxonomic name string.
    This function handles common conventions like 'Candidatus', 'uncultured', 'sp.',
    and also removes square brackets from the name.

    Args:
        name (str): The full taxonomic name (e.g., "[Arthrobacter]", "Microbacterium lemovicicum",
                    "Candidatus Pelagibacter ubique", "uncultured bacterium").

    Returns:
        tuple: A tuple containing (genus, species) strings.
    """
    # Remove square brackets from the name first
    name = name.replace("[", "").replace("]", "").strip()

    parts = name.split()  # Split the name into words

    species = name  # By default, the full name is considered the species
    genus = "Unknown"  # Default genus if not identified

    if not parts:
        # If the name is empty or just whitespace, return defaults
        return genus, species

    # Handle special cases for genus extraction based on common prefixes
    if "Candidatus" in name:
        # For names like "Candidatus Pelagibacter ubique", genus is "Pelagibacter"
        if len(parts) > 1:
            genus = parts[1]
        else:
            genus = "Candidatus"  # If only "Candidatus" is present
    elif "uncultured" in name:
        # For names like "uncultured bacterium" or "uncultured Microbacterium sp."
        # Try to find a more specific genus if available after "uncultured"
        if len(parts) > 1 and parts[1].lower() not in [
            "bacterium",
            "archaeon",
            "organism",
            "clone",
            "fungus",
            "virus",
        ]:
            genus = parts[1]
        elif len(parts) > 2 and parts[2].lower() not in [
            "bacterium",
            "archaeon",
            "organism",
            "clone",
            "fungus",
            "virus",
        ]:
            genus = parts[2]
        else:
            genus = "Uncultured"  # Fallback if no specific genus name follows
    elif "sp." in name and len(parts) > 1:
        # For names like "Microbacterium sp. LWO13-1.2", genus is "Microbacterium"
        genus = parts[0]
    else:
        # For typical binomial names like "Microbacterium lemovicicum", genus is the first word
        genus = parts[0]

    return genus, species


def reindex_filtered_table(df, id_column_name, id_prefix):
    """
    Re-indexes a filtered DataFrame sequentially (e.g., Species1, Species2)
    and converts the ID from index to a regular column.

    Args:
        df (pd.DataFrame): The filtered DataFrame to re-index. Its index
                           is expected to be the original IDs (e.g., 'Species_ID', 'Genus_ID').
        id_column_name (str): The desired name for the re-indexed ID column
                              (e.g., 'Species_ID', 'Genus_ID').
        id_prefix (str): The prefix for the new sequential IDs (e.g., 'Species', 'Genus').

    Returns:
        pd.DataFrame: The re-indexed DataFrame with the ID as a regular column.
                      Returns an empty DataFrame if the input is empty.
    """
    if df.empty:
        return pd.DataFrame()

    # Ensure the DataFrame's index is sorted numerically based on the prefix
    # This is crucial for correct sequential re-indexing (e.g., Species1, Species2, ... not Species1, Species10)
    df = df.sort_index(key=lambda x: x.str[len(id_prefix) :].astype(int))

    # Get the current index (original IDs that passed filtering)
    old_ids = df.index.tolist()

    # Create a new mapping for sequential IDs
    new_id_map = {old_id: f"{id_prefix}{i+1}" for i, old_id in enumerate(old_ids)}

    # Apply the new mapping to the DataFrame's index
    df.rename(index=new_id_map, inplace=True)
    df.index.name = id_column_name  # Ensure index name is maintained

    # Reset the index to convert the ID from index to a regular column
    df.reset_index(inplace=True)

    return df


def process_bracken_data(
    input_dir,
    output_dir,
    metadata_df,
    sample_id_col,
    config,
):
    """
    Processes metagenomics data from Bracken files and generates
    abundance, taxonomy, metadata, and summary tables at both OTU and genus levels.

    Args:
        input_dir (Path): The directory containing the Bracken input files.
        output_dir (Path): The directory where output files for this source will be saved.
        metadata_df (pd.DataFrame): The main metadata DataFrame.
        sample_id_col (str): The name of the column in metadata_df that contains sample IDs.
        config (dict): The loaded configuration dictionary.

    Returns:
        set: A set of sample IDs processed by this source.
    """
    logger.info(f"\n--- Processing Bracken data ---")

    all_species_ids = {}  # Maps taxonomy_id to a unique Species ID
    species_counter = 1
    abundance_data = defaultdict(dict)
    taxonomy_data = []
    processed_samples = set()
    skipped_samples_from_metadata = []

    input_files = list(input_dir.glob(f"*.bracken"))
    logger.info(f"Found {len(input_files)} .bracken files in {input_dir}")

    # Ensure output directory exists for this source
    output_dir.mkdir(parents=True, exist_ok=True)

    # Define keywords to filter out undesirable species names
    undesirable_species_keywords = ["endosymbiont", "secondary", "unidentified"]

    for sample_file in input_files:
        sample_id = sample_file.stem

        # Check if the sample ID exists in the provided metadata
        if sample_id not in metadata_df[sample_id_col].values:
            skipped_samples_from_metadata.append(sample_id)
            continue

        logger.info(f"Processing Bracken file: {sample_file.name}")

        df = pd.read_csv(sample_file, sep="\t")

        # Filter out rows with undesirable species names
        initial_species_count = len(df)
        df = df[
            ~df["name"]
            .astype(str)
            .str.contains("|".join(undesirable_species_keywords), case=False, na=False)
        ]
        filtered_species_count = initial_species_count - len(df)
        if filtered_species_count > 0:
            logger.info(
                f"ℹ️ Filtered out {filtered_species_count} species from sample {sample_id} due to undesirable names."
            )

        for _, row in df.iterrows():
            taxonomy_id = row["taxonomy_id"]
            name = row["name"]
            reads = row["kraken_assigned_reads"]

            if taxonomy_id not in all_species_ids:
                species_id = f"Species{species_counter}"
                all_species_ids[taxonomy_id] = species_id
                species_counter += 1
                genus, species_name = get_genus_and_species(name)
                taxonomy_data.append(
                    {
                        "Species_ID": species_id,
                        "Genus": genus,
                        "Species_Name": species_name,  # Renamed for clarity
                        "taxonomy_id": taxonomy_id,
                    }
                )
            abundance_data[all_species_ids[taxonomy_id]][sample_id] = reads
        processed_samples.add(sample_id)

    logger.info(f"✔ Unique samples processed for Bracken: {len(processed_samples)}")
    logger.info(f"✔ Unique Species detected for Bracken: {len(all_species_ids)}")

    if skipped_samples_from_metadata:
        logger.warning(
            f"⚠ Skipped {len(skipped_samples_from_metadata)} Bracken samples not found in metadata:"
        )
        logger.warning(", ".join(skipped_samples_from_metadata))

    # Create and save Species-level abundance table
    abundance_df = pd.DataFrame.from_dict(abundance_data, orient="index")
    abundance_df.index.name = "Species_ID"  # Renamed index
    abundance_df.fillna(0, inplace=True)
    abundance_df = abundance_df.astype(int)

    if not abundance_df.empty:
        abundance_df.to_csv(output_dir / "species_abundance_table.csv", index=True)
        logger.info(f"📂 Bracken Species-level abundance table saved to {output_dir}.")

        # Filter Species-level abundance
        filtered_abundance_df = filt_otus_by_relfreq_and_prev(
            abundance_df,
            min_rel_freq=config.get("min_rel_freq", 0.00005),
            min_prevalence=config.get("min_prevalence", 0.25),
        )
        # Apply re-indexing for filtered species abundance table
        filtered_abundance_df = reindex_filtered_table(
            filtered_abundance_df, "Species_ID", "Species"
        )

        filtered_abundance_df.to_csv(
            output_dir / "species_abundance_table_filtered.csv", index=False
        )
        logger.info(
            f"📂 Bracken filtered Species-level abundance table saved to {output_dir}."
        )
    else:
        logger.info(
            f"No abundance data generated for Bracken. Skipping abundance table creation."
        )

    # Create and save Species-level taxonomy table
    taxonomy_df = pd.DataFrame(taxonomy_data)
    if not taxonomy_df.empty:
        taxonomy_df.drop_duplicates(
            subset=["Species_ID"], inplace=True
        )  # Use Species_ID for uniqueness
        taxonomy_df.set_index("Species_ID", inplace=True)  # Set Species_ID as index
        taxonomy_df.to_csv(output_dir / "species_taxonomic_table.csv", index=True)
        logger.info(f"📂 Bracken Species-level taxonomy table saved to {output_dir}.")

        # Export filtered Species-level taxonomy table
        if (
            not abundance_df.empty and not filtered_abundance_df.empty
        ):  # Ensure filtered_abundance_df is not empty
            # Filter by Species_ID from the filtered abundance table (which now has Species_ID as a column)
            # IMPORTANT: filtered_abundance_df now has 'Species_ID' as a column, not index.
            # So, filter original taxonomy_df using the 'Species_ID' column of the *re-indexed* filtered_abundance_df
            filtered_taxonomy_df = taxonomy_df[
                taxonomy_df.index.isin(filtered_abundance_df["Species_ID"])
            ]

            # Apply re-indexing for filtered species taxonomy table
            filtered_taxonomy_df = reindex_filtered_table(
                filtered_taxonomy_df, "Species_ID", "Species"
            )

            filtered_taxonomy_df.to_csv(
                output_dir / "species_taxonomic_table_filtered.csv", index=False
            )
            logger.info(
                f"📂 Bracken filtered Species-level taxonomy table saved to {output_dir}."
            )
        else:
            logger.info(
                "No filtered abundance data available to create filtered Species-level taxonomy table."
            )
    else:
        logger.info(
            f"No taxonomy data generated for Bracken. Skipping taxonomy table creation."
        )

    # --- Generate and save Genus-level tables ---
    if not abundance_df.empty and not taxonomy_df.empty:
        # Create Genus-level abundance table
        # Merge abundance with taxonomy to get genus information for aggregation
        # Ensure 'Species_ID' is reset as it's the index of abundance_df
        abundance_with_taxonomy = abundance_df.reset_index().merge(
            taxonomy_df.reset_index()[["Species_ID", "Genus", "taxonomy_id"]],
            on="Species_ID",
            how="left",
        )

        # Create a mapping from Genus name to a unique Genus_ID (e.g., Genus1, Genus2)
        unique_genera = abundance_with_taxonomy["Genus"].dropna().unique()
        # Explicitly sort unique genera to ensure consistent GenusID assignment across runs
        sorted_unique_genera = sorted(unique_genera)
        genus_id_map = {
            genus: str(f"Genus{i+1}") for i, genus in enumerate(sorted_unique_genera)
        }

        # Add Genus_ID column to the merged dataframe
        abundance_with_taxonomy["Genus_ID"] = abundance_with_taxonomy["Genus"].map(
            genus_id_map
        )

        # Group by Genus_ID and sum the read counts for each sample
        # Exclude non-numeric and identifying columns before summing
        cols_to_sum = [
            col
            for col in abundance_with_taxonomy.columns
            if col
            not in ["Species_ID", "Genus", "Genus_ID", "Species_Name", "taxonomy_id"]
        ]

        genus_abundance_df = abundance_with_taxonomy.groupby("Genus_ID")[
            cols_to_sum
        ].sum()
        genus_abundance_df.index.name = "Genus_ID"  # Set Genus_ID as index
        genus_abundance_df = genus_abundance_df.astype(int)

        # Sort the genus abundance DataFrame index numerically for consistent ordering
        genus_abundance_df = genus_abundance_df.sort_index(
            key=lambda x: x.str[len("Genus") :].astype(int)
        )

        # Save unfiltered genus-level abundance table
        genus_abundance_df.to_csv(output_dir / "genus_abundance_table.csv", index=True)
        logger.info(f"📂 Bracken Genus-level abundance table saved to {output_dir}.")

        # Filter genus-level abundance table
        filtered_genus_abundance_df = filt_otus_by_relfreq_and_prev(
            genus_abundance_df,
            min_rel_freq=config.get("min_rel_freq", 0.00005),
            min_prevalence=config.get("min_prevalence", 0.25),
        )
        # Apply re-indexing for filtered genus abundance table
        filtered_genus_abundance_df = reindex_filtered_table(
            filtered_genus_abundance_df, "Genus_ID", "Genus"
        )

        filtered_genus_abundance_df.to_csv(
            output_dir / "genus_abundance_table_filtered.csv", index=False
        )
        logger.info(
            f"📂 Bracken filtered Genus-level abundance table saved to {output_dir}."
        )

        # Create Genus-level taxonomy table with 'Genus_ID', 'Genus_Name', 'taxonomy_id'
        # Get unique genera, their assigned Genus_ID, and their first associated taxonomy_id
        genus_taxonomy_data = []
        for (
            genus_name,
            genus_id,
        ) in genus_id_map.items():  # genus_id_map is already sorted by genus_name
            # Find the first taxonomy_id associated with this genus
            # Ensure we only look for taxonomy_ids from the *original* taxonomy_df that are for this genus
            if (
                genus_name in taxonomy_df["Genus"].values
            ):  # Check if genus_name exists in the original taxonomy_df
                first_taxonomy_id = taxonomy_df[taxonomy_df["Genus"] == genus_name][
                    "taxonomy_id"
                ].iloc[0]
                genus_taxonomy_data.append(
                    {
                        "Genus_ID": str(genus_id),  # Explicitly cast to str
                        "Genus_Name": genus_name,
                        "taxonomy_id": first_taxonomy_id,
                    }
                )

        genus_taxonomy_df = pd.DataFrame(genus_taxonomy_data)
        genus_taxonomy_df.set_index("Genus_ID", inplace=True)  # Set Genus_ID as index

        # Sort the genus taxonomy DataFrame index numerically for consistent ordering
        genus_taxonomy_df = genus_taxonomy_df.sort_index(
            key=lambda x: x.str[len("Genus") :].astype(int)
        )

        genus_taxonomy_df.to_csv(output_dir / "genus_taxonomic_table.csv", index=True)
        logger.info(f"📂 Bracken Genus-level taxonomy table saved to {output_dir}.")

        # Export filtered Genus-level taxonomy table
        if not filtered_genus_abundance_df.empty:
            # Filter genus_taxonomy_df based on the Genus_ID column from filtered_genus_abundance_df
            # IMPORTANT: filtered_genus_abundance_df now has 'Genus_ID' as a column, not index.
            # So, filter original genus_taxonomy_df using the 'Genus_ID' column of the *re-indexed* filtered_genus_abundance_df
            filtered_genus_taxonomy_df = genus_taxonomy_df[
                genus_taxonomy_df.index.isin(filtered_genus_abundance_df["Genus_ID"])
            ]

            # Apply re-indexing for filtered genus taxonomy table
            filtered_genus_taxonomy_df = reindex_filtered_table(
                filtered_genus_taxonomy_df, "Genus_ID", "Genus"
            )

            filtered_genus_taxonomy_df.to_csv(
                output_dir / "genus_taxonomic_table_filtered.csv", index=False
            )
            logger.info(
                f"📂 Bracken filtered Genus-level taxonomy table saved to {output_dir}."
            )
        else:
            logger.info(
                "No filtered genus abundance data available to create filtered Genus-level taxonomy table."
            )
    else:
        logger.info(
            "Skipping Genus-level table creation due to empty abundance or taxonomy data."
        )

    # Create and save metadata table for samples processed by this source
    if processed_samples:
        source_metadata = metadata_df[
            metadata_df[sample_id_col].isin(processed_samples)
        ].copy()
        source_metadata.to_csv(output_dir / "metadata.csv", index=False)
        logger.info(f"📂 Bracken metadata table saved to {output_dir}.")
    else:
        logger.info(
            f"No samples processed for Bracken. Skipping metadata table creation."
        )

    # Create an OTU summary table (using Species_ID for consistency)
    otu_summary = []
    if not abundance_df.empty:
        total_species_in_source = len(all_species_ids)  # Renamed for clarity
        sample_ids_in_abundance = abundance_df.columns.tolist()

        if sample_ids_in_abundance:
            core_species = set(
                abundance_df[abundance_df[sample_ids_in_abundance[0]] > 0].index
            )
        else:
            core_species = set()

        for sample_id in sample_ids_in_abundance:
            sample_species = (abundance_df[sample_id] > 0).sum()  # Renamed for clarity
            percentage = (
                (sample_species / total_species_in_source) * 100
                if total_species_in_source > 0
                else 0
            )

            otu_summary.append(
                {
                    "Sample": sample_id,
                    "Species_in_Sample": sample_species,  # Renamed for clarity
                    "Percentage_of_Source_Species": round(
                        percentage, 2
                    ),  # Renamed for clarity
                    "Total_Species_in_Source": "—",  # Renamed for clarity
                    "Core_Species_in_All_Source_Samples": "—",  # Renamed for clarity
                }
            )
            core_species = core_species.intersection(
                set(abundance_df[abundance_df[sample_id] > 0].index)
            )

        core_species_count = len(core_species)  # Renamed for clarity
        otu_summary.append(
            {
                "Sample": "ALL",
                "Species_in_Sample": "—",
                "Percentage_of_Source_Species": "—",
                "Total_Species_in_Source": total_species_in_source,
                "Core_Species_in_All_Source_Samples": core_species_count,
            }
        )
        summary_df = pd.DataFrame(otu_summary)
        summary_df.to_csv(
            output_dir
            / "otu_summary.csv",  # Filename remains 'otu_summary' as it's a general summary
            index=False,
            encoding="utf-8-sig",
        )
        logger.info(f"📂 Bracken OTU summary saved to {output_dir}.")
    else:
        logger.info(
            f"No abundance data for Bracken OTU summary. Skipping OTU summary creation."
        )

    return processed_samples


def main(config_path):
    """
    Main function to orchestrate the metagenomics data processing.
    It loads configuration and processes Bracken files.

    Args:
        config_path (str): Path to the YAML configuration file.
    """
    global CONFIG, logger

    # Initialize the logger at the very beginning of main
    logger, log_file_path = get_logger(log_suffix="parser_run", display=True)
    logger.info(f"Starting metagenomics data processing. Log file: {log_file_path}")

    CONFIG = load_config(config_path)

    metadata_file = Path(CONFIG["metadata_file"])
    if not metadata_file.exists():
        logger.error(
            f"Error: Metadata file not found at {metadata_file}. Please check your config."
        )
        return

    metadata_df = pd.read_csv(metadata_file)
    sample_id_col = CONFIG["metadata_columns"]["sample_id"]

    # Process Bracken data
    processed_samples = process_bracken_data(
        input_dir=Path(CONFIG["input_dir"]),
        output_dir=Path(CONFIG["output_dir"]),
        metadata_df=metadata_df,
        sample_id_col=sample_id_col,
        config=CONFIG,
    )

    logger.info("\n--- All processing complete ---")
    logger.info(f"Total unique samples processed: {len(processed_samples)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process metagenomics data from Bracken reports."
    )
    parser.add_argument(
        "--config", default="parser_config.yaml", help="Path to config file"
    )
    args = parser.parse_args()

    main(args.config)
