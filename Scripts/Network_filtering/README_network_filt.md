# README: Network Filtering ✂️

This directory contains R scripts for **refining and cleaning microbial interaction networks** after their initial inference. The focus here is on filtering networks based on various thresholds and removing low-connectivity nodes.

---

## Script Descriptions

### Network Cleaning & Filtering

* **`Split_nets_diff_thresholds.R`**: This script **filters an existing `microNet` object's edgelist** across a range of thresholds. It **supports dynamic selection of the filtering column** (`adja`, `asso`, or `diss`). For the `asso` column, it applies a specialized filter that retains edges where the absolute association value is above the threshold (i.e., `value >= threshold` OR `value <= -threshold`). It then **saves new `microNet` objects and corresponding CSV edgelists** for each threshold. This script is intended for single networks.

* **`Split_nets_diff_thresholds_based_on_groups.R`**: This script **filters a combined `microNet` object** (which contains `edgelist1` and `edgelist2` representing different groups). It applies a threshold filter (specifically on the `adja` column by default) to both edgelists and **generates and saves three distinct sets of filtered networks**: a combined network, a PGE-specific network, and a NOPGE-specific network, each with their respective `.rds` and `.csv` files. This script is for networks inferred based on groups.

* **`Remove_singletons_from_net.R`**: This script **loads a network edgelist (CSV)**, converts it into an `igraph` object, and **iteratively removes all "singleton" nodes** (nodes with fewer than 2 connections/edges). It then **saves the cleaned network's edgelist** as a new CSV file. This script is designed to process networks individually.

* **`Remove_singletons_from_biome_nets.R`**: This script **performs the singleton removal process** (removing nodes with fewer than 2 connections) across **multiple biome-specific networks**. It iterates through predefined biomes (`PGE`, `NOPGE`) and various thresholds, cleaning each network identified within the specified folder structure and **saves the cleaned edgelists**. This script is designed to run on multiple network files in a batch.

---

## Running Order for Network Filtering 🚦

Follow these steps to filter and clean your networks:

1.  **Split Networks by Thresholds:**
    * First, run either **`Split_nets_diff_thresholds.R`** (for single networks) or **`Split_nets_diff_thresholds_based_on_groups.R`** (for combined group networks). These scripts will generate multiple versions of your networks, each filtered at different correlation/adjacency thresholds.

2.  **Remove Singletons:**
    * After splitting, run **`Remove_singletons_from_net.R`** (for individual filtered networks) or **`Remove_singletons_from_biome_nets.R`** (for multiple biome-specific filtered networks). These scripts will clean the networks by removing isolated nodes.

---

## General Usage ⚙️

Most R scripts in this folder are executed from the command line using `Rscript`. They utilize the `optparse` package for flexible argument passing.

* For a full list of options and their descriptions for any script using `optparse`, run it with the `--help` flag:

    ```bash
    Rscript your_script_name.R --help
    ```