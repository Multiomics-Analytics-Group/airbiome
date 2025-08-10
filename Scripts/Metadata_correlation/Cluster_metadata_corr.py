#!/usr/bin/env python3
"""
Compute mean CLR-transformed abundances per network community, correlate
module profiles with metadata variables, and generate interactive plots.

This script performs the following steps:
1. Reads abundance, taxonomy, community membership, and metadata files.
2. Applies a Centered Log-Ratio (CLR) transformation to the abundance data.
3. Calculates the mean CLR abundance for each community/module.
4. Correlates these module abundances with specified metadata variables.
5. (Optional) Generates and saves the following interactive plots:
   - Heatmap of mean module abundances.
   - Violin plots of module abundances vs. categorical metadata.
   - Volcano plot of numeric correlations.
   - Heatmap of the top 20 significant numeric correlations.

Usage example:
python module_metadata_corr.py \
  --abundance abundance.csv \
  --taxonomy taxonomy.csv \
  --membership communities.csv \
  --metadata metadata.csv \
  --output-dir results \
  --metadata-cols concentration_pm2_5 Biome \
  --numeric-method spearman \
  --categorical-method kruskal
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Optional, Dict

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr, kruskal, f_oneway
from statsmodels.stats.multitest import multipletests
import plotly.express as px
import plotly.graph_objects as go


def read_csv_file(
    file_path: str, index_col: Optional[int | str] = None
) -> pd.DataFrame:
    """
    Reads a CSV file into a pandas DataFrame.

    Parameters
    ----------
    file_path : str
        Path to the CSV file.
    index_col : int, str or None
        Column to set as index (passed to pd.read_csv). If None, no index column is set.

    Returns
    -------
    pd.DataFrame
        DataFrame with the CSV contents.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    """
    p = Path(file_path)
    if not p.exists():
        raise FileNotFoundError(f"File not found: {file_path}")
    return pd.read_csv(p, index_col=index_col)


def clr_transformation(df: pd.DataFrame, pseudocount: float = 1e-6) -> pd.DataFrame:
    """
    Apply centered log-ratio (CLR) transform to an abundance table.

    Parameters
    ----------
    df : pd.DataFrame
        Abundance table with taxa (species IDs) as rows and samples as columns.
    pseudocount : float
        Small value to add to counts to avoid log(0).

    Returns
    -------
    pd.DataFrame
        CLR-transformed DataFrame (same shape as input).
    """
    arr = df.astype(float) + pseudocount
    log_df = np.log(arr)
    # geometric mean in log-space is mean(log(x))
    gm = log_df.mean(axis=0)
    clr_df = log_df.subtract(gm, axis=1)
    return clr_df


def compute_mean_per_cluster(
    abundance_df: pd.DataFrame,
    taxonomy_df: pd.DataFrame,
    membership_df: pd.DataFrame,
    taxonomy_id_col: Optional[str] = None,
    taxonomy_name_col: Optional[str] = None,
    membership_node_col: str = "Node",
    membership_comm_col: str = "Community",
) -> pd.DataFrame:
    """
    Compute mean CLR abundance for each community (cluster).

    Parameters
    ----------
    abundance_df : pd.DataFrame
        CLR-transformed abundance (rows = species/genus IDs, columns = samples).
    taxonomy_df : pd.DataFrame
        Taxonomy table that contains at least taxonomy_id_col and taxonomy_name_col.
    membership_df : pd.DataFrame
        Node-to-community table that contains membership_node_col and membership_comm_col.
    taxonomy_id_col : str or None
        Column name in taxonomy_df that matches abundance_df.index.
        If None, the function tries to auto-detect a column containing 'id'.
    taxonomy_name_col : str or None
        Column in taxonomy_df that contains taxa names.
        If None, the function tries to auto-detect a column containing 'name'.
    membership_node_col : str
        Column in membership_df with node/species/genus names (default 'Node').
    membership_comm_col : str
        Column in membership_df with community IDs (default 'Community').

    Returns
    -------
    pd.DataFrame
        DataFrame with rows = community IDs and columns = sample IDs (mean CLR abundance).
    """

    # Auto-detect taxonomy_id_col if not provided
    if taxonomy_id_col is None:
        id_cols = [col for col in taxonomy_df.columns if "id" in col.lower()]
        if not id_cols:
            raise ValueError("No ID-like column found in taxonomy dataframe.")
        taxonomy_id_col = id_cols[0]

    # Auto-detect taxonomy_name_col if not provided
    if taxonomy_name_col is None:
        name_cols = [col for col in taxonomy_df.columns if "name" in col.lower()]
        if not name_cols:
            raise ValueError("No Name-like column found in taxonomy dataframe.")
        taxonomy_name_col = name_cols[0]

    # Merge taxonomy with membership on taxonomy_name_col <-> membership_node_col
    merged = taxonomy_df.merge(
        membership_df[[membership_node_col, membership_comm_col]],
        left_on=taxonomy_name_col,
        right_on=membership_node_col,
        how="inner",
    )

    if merged.empty:
        raise ValueError(
            f"No overlap between taxonomy '{taxonomy_name_col}' and membership '{membership_node_col}' columns."
        )

    # Filter species/genus IDs that exist in abundance_df index
    merged = merged[merged[taxonomy_id_col].isin(abundance_df.index)]
    if merged.empty:
        raise ValueError(
            f"After filtering, no IDs from taxonomy '{taxonomy_id_col}' exist in abundance table index."
        )

    # Build community -> taxa list mapping
    comm_to_taxa: Dict[str, List[str]] = (
        merged.groupby(membership_comm_col)[taxonomy_id_col].apply(list).to_dict()
    )

    # Compute mean CLR abundance per community
    rows = {}
    for comm, taxa in comm_to_taxa.items():
        sub = abundance_df.loc[taxa]  # taxa x samples
        rows[comm] = sub.mean(axis=0)  # series length n_samples

    mean_per_cluster = pd.DataFrame(rows).T  # communities x samples
    mean_per_cluster.index.name = str(membership_comm_col)
    return mean_per_cluster


def correlate_modules_metadata(
    mean_per_cluster: pd.DataFrame,
    metadata_df: pd.DataFrame,
    metadata_cols: Optional[List[str]] = None,
    numeric_method: str = "spearman",
    categorical_method: str = "kruskal",
    fdr_per_metadata: bool = True,
) -> pd.DataFrame:
    """
    Correlate each module profile (mean_per_cluster row) with metadata variables.

    Numeric metadata -> Pearson or Spearman (choice via numeric_method).
    Categorical metadata -> Kruskal-Wallis or ANOVA (via categorical_method).

    Parameters
    ----------
    mean_per_cluster : pd.DataFrame
        Rows = communities, columns = sample IDs.
    metadata_df : pd.DataFrame
        Rows = sample IDs (index), columns = metadata variables.
    metadata_cols : list or None
        Optional list of metadata column names to test. If None, all columns are tested.
    numeric_method : str
        'pearson' or 'spearman' (default 'spearman').
    categorical_method : str
        'kruskal' or 'anova' (default 'kruskal').
    fdr_per_metadata : bool
        If True, apply FDR correction separately per metadata variable; otherwise apply global FDR.

    Returns
    -------
    pd.DataFrame
        Results table with columns:
        ['module','metadata','type','stat','p','q','n_samples','n_groups','note']
    """
    # ensure we use only metadata columns requested and exist
    if metadata_cols is None:
        cols_to_test = list(metadata_df.columns)
    else:
        # allow comma-separated single string as well
        if len(metadata_cols) == 1 and "," in metadata_cols[0]:
            cols_to_test = [c.strip() for c in metadata_cols[0].split(",") if c.strip()]
        else:
            cols_to_test = metadata_cols
        missing = [c for c in cols_to_test if c not in metadata_df.columns]
        if missing:
            raise ValueError(
                f"Requested metadata columns not found in metadata file: {missing}"
            )

    # align samples
    samples_common = mean_per_cluster.columns.intersection(metadata_df.index)
    if samples_common.empty:
        raise ValueError(
            "No common sample IDs between module table and metadata index."
        )

    results = []
    skipped_meta = []

    for var in cols_to_test:
        y_raw = metadata_df.loc[samples_common, var]
        # quick constant check: if all identical (after dropping NA) -> skip
        if y_raw.dropna().nunique() <= 1:
            skipped_meta.append((var, "constant"))
            continue

        for module in mean_per_cluster.index:
            x_all = mean_per_cluster.loc[module, samples_common]
            # align x and y and drop NA pairs
            paired = pd.concat([x_all, y_raw], axis=1, join="inner").dropna()
            if paired.shape[0] < 3:
                # not enough paired samples
                note = "too_few_samples"
                results.append(
                    {
                        "module": module,
                        "metadata": var,
                        "type": None,
                        "stat": np.nan,
                        "p": np.nan,
                        "q": np.nan,
                        "n_samples": int(paired.shape[0]),
                        "n_groups": np.nan,
                        "note": note,
                    }
                )
                continue

            x = paired.iloc[:, 0]
            y = paired.iloc[:, 1]

            # Determine if y is numeric (attempt coercion)
            y_numeric = pd.to_numeric(y, errors="coerce")
            n_numeric = y_numeric.notna().sum()
            # if a majority of values can be coerced to numeric and more than 1 unique => treat numeric
            if (
                n_numeric >= max(3, int(0.5 * len(y)))
                and y_numeric.dropna().nunique() > 1
            ):
                # numeric test
                if numeric_method.lower() == "pearson":
                    # check for constant inputs
                    if x.std(ddof=0) == 0 or y_numeric.std(ddof=0) == 0:
                        note = "constant_input"
                        stat_val = np.nan
                        pval = np.nan
                    else:
                        stat_val, pval = pearsonr(x, y_numeric)
                        note = ""
                elif numeric_method.lower() == "spearman":
                    stat_val, pval = spearmanr(x, y_numeric)
                    note = ""
                else:
                    raise ValueError("numeric_method must be 'pearson' or 'spearman'.")

                results.append(
                    {
                        "module": module,
                        "metadata": var,
                        "type": "numeric",
                        "stat": float(stat_val) if np.isfinite(stat_val) else np.nan,
                        "p": float(pval) if np.isfinite(pval) else np.nan,
                        "q": np.nan,  # placeholder, will fill after FDR
                        "n_samples": int(len(x)),
                        "n_groups": np.nan,
                        "note": note,
                    }
                )
            else:
                # categorical test on y
                groups = []
                group_sizes = []
                for g in y.dropna().unique():
                    grp = x[y == g]
                    if grp.size > 0:
                        groups.append(grp.values)
                        group_sizes.append(int(grp.size))
                n_groups = len(groups)
                if n_groups < 2:
                    # cannot test
                    results.append(
                        {
                            "module": module,
                            "metadata": var,
                            "type": "categorical",
                            "stat": np.nan,
                            "p": np.nan,
                            "q": np.nan,
                            "n_samples": int(len(x)),
                            "n_groups": n_groups,
                            "note": "insufficient_groups",
                        }
                    )
                    continue
                try:
                    if categorical_method.lower() == "kruskal":
                        stat_val, pval = kruskal(*groups)
                    elif categorical_method.lower() == "anova":
                        stat_val, pval = f_oneway(*groups)
                    else:
                        raise ValueError(
                            "categorical_method must be 'kruskal' or 'anova'."
                        )
                    note = ""
                except Exception:
                    stat_val, pval = np.nan, np.nan
                    note = "test_error"
                results.append(
                    {
                        "module": module,
                        "metadata": var,
                        "type": "categorical",
                        "stat": float(stat_val) if np.isfinite(stat_val) else np.nan,
                        "p": float(pval) if np.isfinite(pval) else np.nan,
                        "q": np.nan,
                        "n_samples": int(len(x)),
                        "n_groups": n_groups,
                        "note": note,
                    }
                )

    res_df = pd.DataFrame(results)
    if res_df.empty:
        return res_df

    # FDR correction: by metadata variable (default) or global
    if fdr_per_metadata:
        res_df["q"] = np.nan
        for var in res_df["metadata"].unique():
            mask = res_df["metadata"] == var
            pvals = res_df.loc[mask, "p"].fillna(1.0).values
            try:
                qvals = multipletests(pvals, method="fdr_bh")[1]
            except Exception:
                qvals = np.repeat(np.nan, len(pvals))
            res_df.loc[mask, "q"] = qvals
    else:
        pvals = res_df["p"].fillna(1.0).values
        try:
            qvals = multipletests(pvals, method="fdr_bh")[1]
        except Exception:
            qvals = np.repeat(np.nan, len(pvals))
        res_df["q"] = qvals

    # final column order
    cols_out = [
        "module",
        "metadata",
        "type",
        "stat",
        "p",
        "q",
        "n_samples",
        "n_groups",
        "note",
    ]
    res_df = res_df[cols_out]

    return res_df


# --- Plotting Functions ---


def plot_module_means_heatmap(mean_cluster_df: pd.DataFrame, output_dir: Path) -> None:
    """Creates an interactive heatmap of mean module abundances."""
    fig = go.Figure(
        data=go.Heatmap(
            z=mean_cluster_df.values,
            x=mean_cluster_df.columns,
            y=mean_cluster_df.index,
            colorscale="RdBu_r",  # Blue (low) to Red (high)
            hovertemplate="Module: %{y}<br>Sample: %{x}<br>Mean CLR: %{z:.3f}<extra></extra>",
        )
    )
    fig.update_layout(
        template="plotly_white", xaxis_title="Sample ID", yaxis_title="Module"
    )
    fig.update_yaxes(type="category")
    # Save files
    output_path_png = output_dir / "module_means_heatmap.png"
    output_path_json = output_dir / "module_means_heatmap.json"
    fig.write_image(output_path_png, width=1200, height=800)
    fig.write_json(output_path_json, pretty=True)
    print(f"Saved module means heatmap to: {output_path_png} and {output_path_json}")


def plot_module_variability_violins(
    mean_cluster_df: pd.DataFrame, output_dir: Path
) -> None:
    """
    Creates violin plots to show the distribution of CLR values for each module.
    """
    # Transpose to (samples x modules) and melt to a long format
    df_long = mean_cluster_df.T.melt(var_name="Module", value_name="Mean CLR Abundance")

    fig = px.violin(
        df_long,
        x="Module",
        y="Mean CLR Abundance",
        color="Module",
        box=True,
        # Highlight: Change 'points' from False to 'all' to show individual data points
        points="all",
    )
    fig.update_layout(
        template="plotly_white",
        xaxis_title="Module",
        yaxis_title="Mean CLR Abundance",
        showlegend=False,
    )
    fig.update_xaxes(type="category")

    # Save files
    output_path_png = output_dir / "violin_plot_module_variability.png"
    output_path_json = output_dir / "violin_plot_module_variability.json"
    n_modules = len(mean_cluster_df.index)
    fig.write_image(output_path_png, width=max(800, 40 * n_modules), height=600)
    fig.write_json(output_path_json, pretty=True)
    print(
        f"Saved module variability violin plot to: {output_path_png} and {output_path_json}"
    )


def plot_correlation_volcano(
    corr_df: pd.DataFrame, output_dir: Path, q_threshold: float = 0.05
) -> None:
    """Creates an interactive volcano plot for numeric correlation results."""
    numeric_corr = corr_df[corr_df["type"] == "numeric"].copy()
    if numeric_corr.empty:
        print("No numeric correlation results to generate a volcano plot.")
        return

    numeric_corr["-log10(p)"] = -np.log10(numeric_corr["p"].astype(float) + 1e-300)
    numeric_corr["Significance"] = np.where(
        numeric_corr["q"] < q_threshold,
        f"q < {q_threshold}",
        f"q ≥ {q_threshold}",
    )
    numeric_corr["stat"] = pd.to_numeric(numeric_corr["stat"], errors="coerce").fillna(
        0
    )

    fig = px.scatter(
        numeric_corr,
        x="stat",
        y="-log10(p)",
        color="Significance",
        color_discrete_map={
            f"q < {q_threshold}": "#d62728",  # red
            f"q ≥ {q_threshold}": "#7f7f7f",  # grey
        },
        hover_name="module",
        hover_data=["metadata", "p", "q"],
        labels={
            "stat": "Correlation Coefficient (Effect Size)",
            "-log10(p)": "-log₁₀(p-value)",
        },
    )

    fig.add_hline(y=-np.log10(0.05), line_dash="dash", annotation_text="p = 0.05")
    fig.update_layout(template="plotly_white")

    # Save files
    output_path_png = output_dir / "correlation_volcano_plot.png"
    output_path_json = output_dir / "correlation_volcano_plot.json"
    fig.write_image(output_path_png, width=1000, height=800)
    fig.write_json(output_path_json, pretty=True)
    print(f"Saved volcano plot to: {output_path_png} and {output_path_json}")


def plot_correlation_heatmap(
    corr_df: pd.DataFrame, output_dir: Path, n_features: int = 20
) -> None:
    """
    Creates a heatmap of all modules vs. the top N most significant metadata features.
    """
    numeric_corr = corr_df[corr_df["type"] == "numeric"].copy()
    if numeric_corr.empty:
        print("No numeric correlations to generate a heatmap.")
        return
    numeric_corr.dropna(subset=["stat", "q"], inplace=True)

    # Find the top N metadata variables based on their minimum q-value
    min_q_per_meta = numeric_corr.groupby("metadata")["q"].min()
    top_meta_features = min_q_per_meta.nsmallest(n_features).index.tolist()

    if not top_meta_features:
        print(f"No significant metadata features found to create a heatmap.")
        return

    # Filter correlations to only include these top metadata features
    sub_corr_df = numeric_corr[numeric_corr["metadata"].isin(top_meta_features)]

    # Pivot the data: all modules vs. top metadata features
    heatmap_df = sub_corr_df.pivot(
        index="module", columns="metadata", values="stat"
    ).fillna(0)

    # Ensure the columns are ordered by significance
    heatmap_df = heatmap_df[top_meta_features]

    fig = go.Figure(
        data=go.Heatmap(
            z=heatmap_df.values,
            x=heatmap_df.columns,
            y=heatmap_df.index,
            colorscale="RdBu_r",
            zmin=-1,
            zmax=1,
            hovertemplate="Module: %{y}<br>Metadata: %{x}<br>Correlation: %{z:.3f}<extra></extra>",
        )
    )

    fig.update_layout(
        template="plotly_white",
        xaxis_title="Top Metadata Features (by significance)",
        yaxis_title="Module",
    )
    # Explicitly set axis types to handle numeric-like labels correctly
    fig.update_xaxes(type="category")
    fig.update_yaxes(type="category")

    output_path_png = output_dir / "top_correlations_heatmap.png"
    output_path_json = output_dir / "top_correlations_heatmap.json"
    fig.write_image(
        output_path_png, width=1200, height=max(600, 20 * len(heatmap_df.index))
    )
    fig.write_json(output_path_json, pretty=True)
    print(
        f"Saved top correlations heatmap to: {output_path_png} and {output_path_json}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Module <-> metadata correlation pipeline with plotting"
    )
    parser.add_argument(
        "--abundance",
        required=True,
        help="Abundance CSV (rows = species_id, cols = samples)",
    )
    parser.add_argument(
        "--taxonomy",
        required=True,
        help="Taxonomy CSV (must contain Species_ID and Species_Name)",
    )
    parser.add_argument(
        "--membership",
        required=True,
        help="Node-to-community CSV (must contain Node and Community)",
    )
    parser.add_argument(
        "--metadata",
        required=True,
        help="Metadata CSV (rows = samples; index column selectable)",
    )
    parser.add_argument(
        "--output-dir", required=True, help="Directory where outputs will be saved"
    )
    parser.add_argument(
        "--metadata-cols",
        nargs="*",
        default=None,
        help="Optional list of metadata columns to test (space-separated) or comma-separated single string",
    )
    parser.add_argument(
        "--abundance-index",
        type=int,
        default=0,
        help="Index column for abundance CSV (default 0)",
    )
    parser.add_argument(
        "--metadata-index",
        type=int,
        default=0,
        help="Index column for metadata CSV (default 0)",
    )
    parser.add_argument(
        "--numeric-method",
        choices=["pearson", "spearman"],
        default="spearman",
        help="Method for numeric metadata correlation (default: spearman)",
    )
    parser.add_argument(
        "--categorical-method",
        choices=["kruskal", "anova"],
        default="kruskal",
        help="Method for categorical metadata (default: kruskal)",
    )
    parser.add_argument(
        "--fdr-per-metadata",
        action="store_true",
        default=True,
        help="Apply FDR correction per metadata variable (default).",
    )
    parser.add_argument(
        "--skip-plots",
        action="store_true",
        help="Skip generating and saving plots.",
    )
    args = parser.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Read files
    abundance_df = read_csv_file(args.abundance, index_col=args.abundance_index)
    taxonomy_df = read_csv_file(args.taxonomy, index_col=None)
    membership_df = read_csv_file(args.membership, index_col=None)
    metadata_df = read_csv_file(args.metadata, index_col=args.metadata_index)

    # CLR transform
    clr_df = clr_transformation(abundance_df)

    # Compute mean per cluster
    mean_cluster = compute_mean_per_cluster(
        clr_df,
        taxonomy_df,
        membership_df,
    )

    # Save module means (communities x samples)
    module_means_path = outdir / "module_means.csv"
    mean_cluster.to_csv(module_means_path)
    print(f"Saved module means to: {module_means_path}")

    # Correlate with metadata
    res_df = correlate_modules_metadata(
        mean_cluster,
        metadata_df,
        metadata_cols=args.metadata_cols,
        numeric_method=args.numeric_method,
        categorical_method=args.categorical_method,
        fdr_per_metadata=args.fdr_per_metadata,
    )

    corr_path = outdir / "module_metadata_correlations.csv"
    res_df.to_csv(corr_path, index=False)
    print(f"Saved correlation results to: {corr_path}")

    # Write skipped metadata (constant) to a small log file
    if args.metadata_cols is None:
        requested_cols = list(metadata_df.columns)
    elif len(args.metadata_cols) == 1 and "," in args.metadata_cols[0]:
        requested_cols = [
            c.strip() for c in args.metadata_cols[0].split(",") if c.strip()
        ]
    else:
        requested_cols = list(args.metadata_cols)

    skipped = []
    for var in requested_cols:
        if var not in metadata_df.columns:
            skipped.append((var, "missing"))
        elif metadata_df[var].dropna().nunique() <= 1:
            skipped.append((var, "constant"))

    if skipped:
        skipped_path = outdir / "skipped_metadata.txt"
        with open(skipped_path, "w") as fh:
            fh.write("metadata\treason\n")
            for v, reason in skipped:
                fh.write(f"{v}\t{reason}\n")
        print(f"Some metadata variables were skipped. See {skipped_path}")

    # Print small summary and export it as a CSV
    if not res_df.empty:
        sig = res_df[res_df["q"].le(0.05)].sort_values("q").head(20)
        if not sig.empty:
            print("\nTop significant module <-> metadata associations (q <= 0.05):")
            print(
                sig[["module", "metadata", "type", "stat", "p", "q"]].to_string(
                    index=False
                )
            )
        else:
            print("\nNo associations with q <= 0.05 found.")
    else:
        print("\nNo correlation tests were performed.")

    sig_corr_df = res_df[res_df["q"].le(0.05)]
    if not sig_corr_df.empty:
        sig_corr_path = outdir / "significant_module_metadata_corr.csv"
        sig_corr_df.to_csv(sig_corr_path, index=False)
        print(f"Saved significant correlations (q <= 0.05) to: {sig_corr_path}")

    # --- Plotting ---
    if not args.skip_plots:
        print("\n--- Generating Plots ---")
        print(
            "NOTE: Plot export to PNG requires the 'kaleido' package (pip install kaleido)."
        )

        # Plot 1: Heatmap of all module means
        plot_module_means_heatmap(mean_cluster, outdir)

        # Plot 2: Violin plots showing variability of each module
        plot_module_variability_violins(mean_cluster, outdir)

        # Plot 3: Volcano plot for numeric correlations
        plot_correlation_volcano(res_df, outdir)

        # Plot 4: Heatmap for top 20 numeric correlations
        plot_correlation_heatmap(res_df, outdir, n_features=10)


if __name__ == "__main__":
    main()
