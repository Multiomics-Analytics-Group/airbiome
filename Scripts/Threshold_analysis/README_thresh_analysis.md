# README: Threshold Analysis 📊

This directory contains Python scripts for performing **threshold analysis of microbial association networks**, including the **calculation of topological metrics** and **visualization of these metrics** across different thresholds. The scripts are provided in pairs: one set for single datasets and another for processing multiple biomes.

---

## Script Descriptions 💡

### Network Metric Calculation

* **`Calculate_network_metrics_and_get_best_net.py`**:

    * **Purpose**: Computes and consolidates network topological metrics (e.g., nodes, edges, modularity) for networks originating from a **single base directory's threshold-filtered outputs**. It identifies a "best" network based on specific criteria.

    * **Outputs**: A CSV file containing all calculated metrics and a separate CSV file with node clustering data for the identified "best" network.

* **`Calculate_network_metrics_and_get_best_net_bybiome.py`**:

    * **Purpose**: Processes network metrics across **multiple biome subdirectories**, performing the same calculations as its single-biome counterpart. It then consolidates the metrics from all biomes into a single unified CSV.

    * **Outputs**: A unified metrics CSV for all biomes, and individual clustering data CSVs for each biome's "best" network.

### Plotting Scripts

* **`Plot_threshold_analysis_metrics.py`**:

    * **Purpose**: Generates **line plots** visualizing how key network metrics (e.g., Nodes, Edges, Communities, Modularity, Average Degree, Average Clustering Coefficient) change with the association threshold. This script is intended for metrics from a **single dataset**.

    * **Outputs**: PNG and SVG plots.

* **`Plot_threshold_analysis_metrics_bybiome.py`**:

    * **Purpose**: Generates **comparative line plots** from a unified metrics CSV that contains data from **multiple biomes**. Each biome is represented by a distinct line, allowing for easy comparison of trends across thresholds.

    * **Outputs**: PNG and SVG plots.

---

## How to Run the Scripts 🚀

All scripts are executed via the command line using `python` and `argparse` arguments.

### Running Metric Calculation

* **For a single biome/dataset**:

    ```bash
    python Calculate_network_metrics_and_get_best_net.py \
      --base_dir "path/to/network_inference/single_biome_output" \
      --output_dir "path/to/analysis_output/single_biome_results" \
      --network_type "cclasso" \
      --weight_column "adja"
    ```

* **For multiple biomes**:

    ```bash
    python Calculate_network_metrics_and_get_best_net_bybiome.py \
      --base_dir "path/to/network_inference/all_biomes_output" \
      --output_dir "path/to/analysis_output/all_biomes_results" \
      --network_type "cclasso" \
      --weight_column "asso"
    ```

### Running Plotting Scripts

* **For a single biome/dataset's metrics**:

    ```bash
    python Plot_threshold_analysis_metrics.py \
      --input_file "path/to/metrics_csv/single_biome_metrics.csv" \
      --output_dir "path/to/plots/single_biome" \
      --dataset_name "My_Single_Biome_Name"
    ```

* **For multiple biomes' metrics**:

    ```bash
    python Plot_threshold_analysis_metrics_bybiome.py \
      --input_file "path/to/metrics_csv/all_biomes_metrics.csv" \
      --output_dir "path/to/plots/all_biomes" \
      --dataset_name "All_Biomes_Combined"
    ```

Remember to replace `"path/to/..."` with your actual file system paths.