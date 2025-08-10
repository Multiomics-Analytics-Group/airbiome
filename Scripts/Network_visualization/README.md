# README: Network Visualization 🕸️

This directory contains R scripts for **visualizing microbial interaction networks in Cytoscape**. These scripts facilitate the import and mapping of network properties, including node degree, community clusters, and edge weights (with handles signed associations).

---

## Script Descriptions 💡

### 1. `Create_cytoscape_viz.R`

* **Purpose**: This script is designed to import and style a **single network** (from a `.csv` edgelist) into Cytoscape. It applies visual mappings for node size (based on degree), node color (by cluster/community), and edge color/width (based on the specified weight column, handling signed 'asso' values appropriately).
* **Inputs**:
    * A single network `.csv` file (must be "nosinglt" filtered).
    * The desired edge weight column (`adja`, `asso`, or `diss`).
    * A name for the Cytoscape collection.
* **Outputs**: A network within your running Cytoscape instance.

### 2. `Create_cytoscape_viz_bybiome.R`

* **Purpose**: This script automates the visualization process for **multiple networks** across various biomes and thresholds. It loops through specified biome directories and threshold folders, importing and styling each network individually into Cytoscape. This is ideal for batch processing and comparing networks across different conditions.
* **Inputs**:
    * A root directory containing biome-specific subfolders, which in turn contain network files.
    * Network type (e.g., `cclasso`, `spring`).
    * Lists of biomes and thresholds to process.
    * The desired edge weight column (`adja`, `asso`, or `diss`).
    * A prefix for Cytoscape collection names (each biome gets its own collection).
* **Outputs**: Multiple networks organized by biome and threshold within your running Cytoscape instance.

---

## How to Run the Scripts 🚀 

Ensure **Cytoscape Desktop is running** before executing these R scripts.

### 1. Running `Create_cytoscape_viz.R` (Single Network)

To visualize a single network:

```bash
Rscript Create_cytoscape_viz.R \
  --net_path "/path/to/your/network/file_nosinglt.csv" \
  --weight_column "adja" \
  --collection_name "My_Single_Network_Collection"
```

* Replace `/path/to/your/network/file_nosinglt.csv` with the **full path** to your network file.
* Adjust `--weight_column` (`adja`, `asso`, or `diss`) as needed.
* Customize `--collection_name` for organization in Cytoscape.

### 2. Running `Create_cytoscape_viz_bybiome.R` (Multiple Networks/Biomes)

To visualize multiple networks across biomes and thresholds:

```bash
Rscript Create_cytoscape_viz_bybiome.R \
  --root_path "/path/to/your/indiv_biomes_root_folder" \
  --network_type "cclasso" \
  --biomes "BiomeA,BiomeB" \
  --thresholds "0.1:0.7:0.05" \
  --weight_column "asso" \
  --collection_name_prefix "Project_Batch_Viz"
```

* Replace `/path/to/your/indiv_biomes_root_folder` with the top-level directory containing your biome subfolders.
* `--network_type`: Specify the network inference method (e.g., `cclasso`, `spring`).
* `--biomes`: Provide a comma-separated list of biome folder names.
* `--thresholds`: Can be a comma-separated list (e.g., `"0.1,0.2,0.3"`) or a range (`"start:end:step"`, e.g., `"0.1:0.7:0.05"`).
* `--weight_column`: Choose `adja`, `asso`, or `diss`.
* `--collection_name_prefix`: Networks will be organized into collections like `YOUR_PREFIX_BIOMENAME`.

Remember to have the `RCy3`, `igraph`, `dplyr`, `optparse`, and `paletteer` R packages installed.