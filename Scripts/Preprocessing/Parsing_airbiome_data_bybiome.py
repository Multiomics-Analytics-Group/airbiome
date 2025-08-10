import pandas as pd
from collections import defaultdict
from pathlib import Path
import yaml
import argparse


def load_config(config_path):
    """Load configuration from YAML file without path conversion"""
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def filt_otus_by_relfreq_and_prev(
    abundance_df, min_rel_freq=0.00005, min_prevalence=0.25
):
    """
    Filter OTUs by relative frequency and prevalence thresholds.
    """
    rel_abundance = abundance_df.div(abundance_df.sum(axis=0), axis=1)
    presence_matrix = rel_abundance >= min_rel_freq
    prevalence = presence_matrix.sum(axis=1) / presence_matrix.shape[1]
    keep_otus = prevalence[prevalence >= min_prevalence].index
    return abundance_df.loc[keep_otus]


def process_biome(biome_name, biome_config, metadata_df):
    """Process samples for a single biome with late path conversion"""
    # Convert to Path objects exactly when needed
    input_dir = Path(biome_config["input_dir"])
    output_dir = Path(biome_config["output_dir"])

    # Create output directory if needed
    output_dir.mkdir(parents=True, exist_ok=True)

    # Initialize data structures
    all_otus = {}  # {taxonomy_id: otu_id}
    otu_counter = 1
    abundance_data = defaultdict(dict)
    taxonomy_data = []
    sample_species_counts = {}
    skipped_samples = []
    otu_summary = []

    # Get sample list for this biome
    biome_samples = set(
        metadata_df[metadata_df[CONFIG["metadata_columns"]["biome"]] == biome_name][
            CONFIG["metadata_columns"]["sample_id"]
        ].tolist()
    )

    # Process each sample file
    for sample_file in input_dir.glob("*.bracken"):
        sample_id = sample_file.stem

        if sample_id not in biome_samples:
            skipped_samples.append(sample_id)
            continue

        df = pd.read_csv(sample_file, sep="\t")
        sample_species_counts[sample_id] = len(df)

        for _, row in df.iterrows():
            taxonomy_id = row["taxonomy_id"]

            if taxonomy_id not in all_otus:
                otu_id = f"OTU{otu_counter}"
                all_otus[taxonomy_id] = otu_id
                otu_counter += 1

                name = row["name"].strip()
                parts = name.split()

                # Handle special cases for genus extraction
                if "Candidatus" in name or "uncultured" in name:
                    genus = parts[1] if len(parts) > 1 else "Unknown"
                else:
                    genus = parts[0] if len(parts) > 0 else "Unknown"

                taxonomy_data.append(
                    {
                        "OTU": otu_id,
                        "Genus": genus if genus else "Unknown",
                        "Species": name,
                        "taxonomy_id": taxonomy_id,
                    }
                )

            abundance_data[all_otus[taxonomy_id]][sample_id] = row[
                "kraken_assigned_reads"
            ]

    # Validation summary
    print(f"✔ Processed {len(biome_samples)} samples for biome: {biome_name}")
    print(f"✔ Unique OTUs detected: {len(all_otus)}")

    if skipped_samples:
        print(f"⚠ Skipped {len(skipped_samples)} samples not found in metadata:")
        print(", ".join(skipped_samples))

    # Check species counts match
    for sample_id, expected_species in sample_species_counts.items():
        actual_species = sum(sample_id in v for v in abundance_data.values())
        if actual_species != expected_species:
            print(
                f"⚠ Mismatch in species count for sample {sample_id}: "
                f"{actual_species} added vs {expected_species} in file"
            )

    # Create and save abundance table
    abundance_df = pd.DataFrame.from_dict(abundance_data, orient="index")
    abundance_df.index.name = "OTU"
    abundance_df.fillna(0, inplace=True)
    abundance_df = abundance_df.astype(int)

    # Save unfiltered table (optional)
    abundance_df.to_csv(output_dir / f"{biome_name}_abundance_table.csv", index=True)
    print(f"📂 Abundance table saved for {biome_name}.")

    # Filter OTUs based on relative frequency and prevalence
    filtered_abundance_df = filt_otus_by_relfreq_and_prev(
        abundance_df,
        min_rel_freq=CONFIG.get("min_rel_freq", 0.00005),
        min_prevalence=CONFIG.get("min_prevalence", 0.25),
    )

    # Save filtered abundance table
    filtered_abundance_df.to_csv(
        output_dir / f"{biome_name}_abundance_table_filtered.csv", index=True
    )
    print(f"📂 Filtered abundance table saved for {biome_name}.")

    # Create and save taxonomy table
    taxonomy_df = pd.DataFrame(taxonomy_data)
    taxonomy_df.set_index("OTU", inplace=True)
    taxonomy_df.to_csv(output_dir / f"{biome_name}_taxonomic_table.csv", index=True)
    print(f"📂 Taxonomy table saved for {biome_name}.")

    # Create and save metadata table
    biome_metadata = metadata_df[metadata_df["Library"].isin(biome_samples)].copy()
    biome_metadata.to_csv(output_dir / f"{biome_name}_metadata.csv", index=False)
    print(f"📂 Metadata table saved for {biome_name}.")

    # Create an OTU summary table
    total_otus_in_biome = len(all_otus)
    sample_ids = abundance_df.columns.tolist()
    core_otus = set(abundance_df[abundance_df[sample_ids[0]] > 0].index)

    # Iterate through samples to gather data
    for sample_id in sample_ids:
        # Count OTUs for this sample
        sample_otus = (abundance_df[sample_id] > 0).sum()
        percentage = (sample_otus / total_otus_in_biome) * 100

        otu_summary.append(
            {
                "Sample": sample_id,
                "OTUs_in_Sample": sample_otus,
                "Percentage_of_Biome_OTUs": round(percentage, 2),
                "Total_OTUs_in_Biome": "—",
                "Core_OTUs_in_All_Samples": "—",
            }
        )

        # Update core OTUs (only OTUs present in this sample)
        core_otus = core_otus.intersection(
            set(abundance_df[abundance_df[sample_id] > 0].index)
        )

    # Calculate Core OTUs present in all samples
    core_otus_count = len(core_otus)

    # Add the "ALL" sample summary
    otu_summary.append(
        {
            "Sample": "ALL",
            "OTUs_in_Sample": "—",
            "Percentage_of_Biome_OTUs": "—",
            "Total_OTUs_in_Biome": total_otus_in_biome,
            "Core_OTUs_in_All_Samples": core_otus_count,
        }
    )

    # Write the summary to a CSV file
    summary_df = pd.DataFrame(otu_summary)
    summary_df.to_csv(
        output_dir / f"{biome_name}_otu_summary.csv",
        index=False,
        encoding="utf-8-sig",
    )
    print(f"📂 OTU summary saved for {biome_name}.")


def main(config_path):
    global CONFIG
    CONFIG = load_config(config_path)
    metadata_df = pd.read_csv(CONFIG["metadata_file"])

    for biome_name, biome_config in CONFIG["biomes"].items():
        print(f"Processing {biome_name}...")
        process_biome(biome_name, biome_config, metadata_df)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config", default="parser_config.yaml", help="Path to config file"
    )
    args = parser.parse_args()

    main(args.config)
