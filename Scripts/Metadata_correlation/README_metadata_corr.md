# README: Module-Metadata Correlation 📊

This directory contains the `Cluster_metadata_corr.py` script, a Python tool for analyzing the relationship between microbial community modules and metadata variables. 

-----

## Script Description 💡

### `Cluster_metadata_corr.py`

  * **Purpose**: This script is designed to automate the process of correlating **mean module abundances** with metadata. It handles data preprocessing (CLR transformation), performs statistical tests (Spearman, Pearson, Kruskal-Wallis, ANOVA) on both numeric and categorical variables, applies FDR correction, and outputs comprehensive results and plots.

  * **Inputs**:
    * `--abundance`: A CSV file of microbial abundances (Taxon IDs as rows, sample IDs as columns).
    * `--taxonomy`: A CSV file linking taxon IDs to names (The script can auto-detect "id" and "name" columns).
    * `--membership`: A CSV file mapping taxa to community modules.
    * `--metadata`: A CSV file containing sample metadata (sample IDs as the index and metadata variables as columns.).

  * **Outputs**:
    * **`module_means.csv`**: Mean CLR abundance per module.
    * **`module_metadata_correlations.csv`**: A table of all correlation results.
    * **`significant_module_metadata_corr.csv`**: Filtered results with a q-value ≤ 0.05.
    * **`skipped_metadata.txt`**: A log file of any skipped metadata columns.
    * **Plots**: Interactive `.json` and static `.png` versions of all generated plots.

-----

## Features ✨

  * **Data Transformation**: Applies a **Centered Log-Ratio (CLR)** transformation to abundance data to handle its compositional nature.
  * **Module Profiling**: Computes the mean CLR abundance for each community module across all samples.
  * **Statistical Analysis**:
      * **Numeric Metadata**: Uses **Spearman** (default) or **Pearson** correlation.
      * **Categorical Metadata**: Uses **Kruskal-Wallis** (default) or **ANOVA**.
      * Includes **FDR correction** for multiple comparisons.
  * **Visualization**: Generates interactive plots (`.json` and `.png`) for deep data exploration.

-----

## Requirements 🛠️

Install the necessary Python packages, including `kaleido` for exporting static plots:

```bash
pip install pandas numpy scipy statsmodels plotly kaleido
```

-----

## How to Run the Script 🚀

To analyze module-metadata correlations and generate plots:

```bash
python module_metadata_corr.py \
  --abundance /path/to/your/abundance.csv \
  --taxonomy /path/to/your/taxonomy.csv \
  --membership /path/to/your/communities.csv \
  --metadata /path/to/your/metadata.csv \
  --output-dir results \
  --metadata-cols concentration_pm2_5 Biome \
  --numeric-method spearman \
  --categorical-method kruskal
```

  * Replace the `/path/to/your/...` placeholders with the **paths** to your data files.
  * `--metadata-cols`: Provide a space-separated list of metadata columns to test. Omit this argument to test all columns.
  * `--numeric-method`: Choose `spearman` or `pearson` for numeric correlations.
  * `--categorical-method`: Choose `kruskal` or `anova` for categorical comparisons.

Remember to have `pandas`, `numpy`, `scipy`, `statsmodels`, `plotly`, and `kaleido` installed.
