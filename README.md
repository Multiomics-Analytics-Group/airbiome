# Microbial Network Analysis Workflow 🦠📊

This repository contains R and Python scripts organized into modules to perform a comprehensive microbial association network (MAN) analysis, 
from raw data preprocessing to network comparisons and metadata correlations. Here, we preprocess abundance and metadata files, build MANs
using the [NetComi R package](https://github.com/stefpeschel/NetCoMi/), filter and analyze these networks, visualize them in Cytoscape, 
and correlate network modules with metadata variables.

---

## 🚀 Workflow Overview

The general workflow for analyzing MANs in this repository follows these steps:

1.  **Preprocessing**: Prepare raw abundance and metadata for analysis. 

2.  **Network Inference**: Construct microbial association networks.

3.  **Network Filtering**: Split networks based on various thresholds and remove singletons or low-connectivity nodes.

4.  **Threshold Analysis**: Evaluate network topological metrics across networks with thresholds and identify optimal networks 
based on modularity and average clsutering coefficient.

5.  **Network Visualization**: Generate interactive visualizations of networks in Cytoscape.

6.  **Network Analysis**: Perform network analysis using NetComi and custom Python scripts.

7.  **Metadata Correlation**: Investigate relationships between network modules and metadata variables.

8.  **Networks Comparison**: Compare network structures, overlaps, and clustering similarities.

## 📁 Scripts & Their Purpose

Each subfolder within the `Scripts` directory contains specialized tools. For detailed usage instructions and parameters for individual scripts, refer to the `README.md` file located within each respective subfolder.

| Folder | Purpose | Key Scripts | Outputs |
| :----------------------- | :------------------------------------------------------------------------------------------------------- | :-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | :--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `Preprocessing` | Initial data handling: parsing raw data into usable formats and exploratory data analysis (EDA). | `Parsing_airbiome_data_allsamples.py`, `Parsing_airbiome_data_bybiome.py`, `EDA_all_samples.ipynb`, `EDA_bybiome.ipynb` | Cleaned abundance, taxonomic, and metadata files, EDA plots. |
| `Network_inference` | Building microbial interaction networks using various methods. | `Build_microb_net.R`, `Build_microb_net_sp_gen.R`, `Build_microb_nets_based_on_groups.R`, `Create_phyloseq_object.R`, `Create_phyloseq_objects_bybiome.R` | `microNet` R objects and corresponding edgelist CSVs, often saved by network type and threshold. |
| `Threshold_analysis` | Calculating and visualizing topological metrics of networks across different thresholds. | `Calculate_network_metrics_and_get_best_net.py` (single dataset), `Calculate_network_metrics_and_get_best_net_bybiome.py` (multiple biomes), `Plot_threshold_analysis_metrics.py` (plotting single dataset), `Plot_threshold_analysis_metrics_bybiome.py` (plotting multiple biomes). | CSV files with calculated network metrics, PNG/SVG plots illustrating metric trends, and node clustering data for "best" networks. |
| `Network_filtering` | Refining and cleaning inferred networks by applying thresholds and removing low-connectivity nodes. | `Split_nets_diff_thresholds.R`, `Split_nets_diff_thresholds_based_ongroups.R`, `Remove_singletons_from_net.R`, `Remove_singletons_from_biome_nets.R` | Filtered `microNet` objects and cleaned edgelist CSVs (typically with `nosinglt` in the filename). |
| `Network_visualization` | Creating dynamic visualizations of networks in Cytoscape. | `Create_cytoscape_viz.R` (single network), `Create_cytoscape_viz_bybiome.R` (multiple networks/biomes). | Networks imported and visually styled within a running Cytoscape instance, mapping node properties (e.g., degree, cluster) and edge properties (e.g., weight, sign). |
| `Network_analysis` | Advanced network analyses. | `Centrality_analysis.py`, `Multi_layer_network.py`, `Network_analysis_netcomi.R`, `Network_analysis_netcomi_bybiome.R` | Various analytical results, potentially including centrality measures, multi-layer network structures, and other network-specific metrics. |
| `Metadata_correlation` | Correlating microbial community modules with sample metadata. | `Cluster_metadata_corr.py` | CSV reports of module-metadata correlations (including statistical tests and FDR correction), and corresponding plots. |
| `Networks_comparison` | Comparing different networks. | `Abund_analysis_overlap_and_unique_sp.ipynb`, `Calculate_clustering_similarity_2nets.py`, `Find_overlap_and_unique_otus_2nets.py`, `Find_overlapping_edges_2nets.py` | Metrics and reports on network overlap, unique features, and clustering similarities. |

For detailed instructions on running each script and understanding specific outputs, please refer to the individual `README.md` files located within each script subfolder (e.g., `Scripts/Network_filtering/README.md`).

## 📦 Output Folders

The `Output` directory is structured to store results from various analysis stages:

* **`Best_nets`**: Stores information about optimal networks identified during threshold analysis.

* **`EDA`**: Contains results from initial exploratory data analysis.

* **`Network_inference`**: Holds the raw and intermediate results of network inference.

* **`Threshold_analysis`**: Stores calculated metrics and plots from threshold analysis.

## Summary Reports

We created sumamry reports with the main results of this project using [VueGen](https://github.com/Multiomics-Analytics-Group/vuegen), a tool to automate the generation of scientific reports:  

[![HTML5](https://img.shields.io/badge/html5-%23E34F26.svg?style=for-the-badge&logo=html5&logoColor=white)][html-report]
[![Streamlit](https://img.shields.io/badge/Streamlit-%23FE4B4B.svg?style=for-the-badge&logo=streamlit&logoColor=white)][streamlit-report]   

Further details on the reports and the source code to generate them is available on this [GitHub repository][report-github-repo].

[html-report]: https://multiomics-analytics-group.github.io/airbiome_microb_asso_net_report/
[streamlit-report]: https://airbiome-microb-ass-net-summary.streamlit.app/
[report-github-repo]: https://github.com/Multiomics-Analytics-Group/airbiome_microb_asso_net_report
