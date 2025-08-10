import pandas as pd
import networkx as nx
import numpy as np
import os
import argparse
from typing import Dict, Any, List, Tuple, Union


def create_networkx_graph(file_path: str, weight_column: str) -> nx.Graph:
    """
    Reads an edge list from a CSV file and creates a NetworkX graph.

    Parameters
    ----------
    file_path : str
        The path to the CSV file containing the edge list.
    weight_column : str
        The name of the column in the CSV file that should be used as edge weight.

    Returns
    -------
    nx.Graph
        A NetworkX graph object created from the edge list.

    Raises
    ------
    ValueError
        If the specified `weight_column` does not exist in the CSV file.
    """
    df = pd.read_csv(file_path)

    # Validate if the weight_column exists in the DataFrame
    if weight_column not in df.columns:
        raise ValueError(
            f"Weight column '{weight_column}' not found in the CSV file: {file_path}"
        )

    G = nx.from_pandas_edgelist(
        df, source="v1", target="v2", edge_attr=weight_column, create_using=nx.Graph()
    )
    return G


def transform_weights_to_absolute(
    G: nx.Graph, weight_attribute: str = "asso"
) -> nx.Graph:
    """
    Transforms the weights of a NetworkX graph to their absolute values for a specified edge attribute.

    This is particularly useful for metrics that do not support signed networks (e.g., modularity,
    clustering coefficient) when a column that can have negative values is used as weight.

    Parameters
    ----------
    G : nx.Graph
        The input NetworkX graph.
    weight_attribute : str, optional
        The name of the edge attribute whose values should be converted to absolute.
        Defaults to 'asso'.

    Returns
    -------
    nx.Graph
        A new NetworkX graph object with the specified edge attribute values transformed to absolute values.
    """
    G_abs = G.copy()
    for u, v, d in G_abs.edges(data=True):
        if weight_attribute in d:
            d[weight_attribute] = abs(d[weight_attribute])
    return G_abs


def calculate_network_metrics(
    G: nx.Graph, weight_column: str
) -> Dict[str, Union[int, float, str]]:
    """
    Calculates general topological metrics for a given NetworkX graph.

    Parameters
    ----------
    G : nx.Graph
        The NetworkX graph object for which to calculate metrics.
    weight_column : str
        The name of the edge attribute to use as weight for calculations (e.g., modularity, average degree).

    Returns
    -------
    Dict[str, Union[int, float, str]]
        A dictionary containing various network topological metrics:
        - Nodes (int): Number of nodes.
        - Edges (int): Number of edges.
        - Number_communities (int or str): Number of communities detected, or 'N/A' if detection fails.
        - Modularity (float or str): Modularity score, or 'N/A' if calculation fails.
        - Avg_degree (float): Average node degree.
        - Avg_clustering_coefficient (float): Average clustering coefficient.
        - Density (float): Graph density.
    """
    num_nodes: int = G.number_of_nodes()
    num_edges: int = G.number_of_edges()

    # Community detection using Louvain method
    if num_nodes == 0:
        return {
            "Nodes": 0,
            "Edges": 0,
            "Number_communities": 0,
            "Modularity": 0,
            "Avg_degree": 0,
            "Avg_clustering_coefficient": 0,
            "Density": 0,
        }

    num_communities: Union[int, str] = "N/A"
    modularity: Union[float, str] = "N/A"
    try:
        communities = list(
            nx.algorithms.community.louvain_communities(
                G, weight=weight_column, seed=12
            )
        )
        num_communities = len(communities)
        modularity = nx.algorithms.community.modularity(
            G, communities, weight=weight_column
        )
    except Exception as e:
        print(f"Warning: Community detection or modularity calculation failed: {e}")

    # Calculate average degree and density
    degrees = [deg for node, deg in G.degree(weight=weight_column)]
    avg_degree = np.mean(degrees)
    density = nx.density(G)

    # Calculate the average clustering coefficient
    avg_clustering_coefficient: float = nx.average_clustering(G, weight=weight_column)

    return {
        "Nodes": num_nodes,
        "Edges": num_edges,
        "Number_communities": num_communities,
        "Modularity": modularity,
        "Avg_degree": avg_degree,
        "Avg_clustering_coefficient": avg_clustering_coefficient,
        "Density": density,
    }


def create_metrics_inconsistent(
    threshold: Union[str, float],
    network_type: str,
    dataset_name: str,
    message: str,
    n_nodes: Union[int, str, None] = None,
    n_edges: Union[int, str, None] = None,
) -> Dict[str, Union[str, int, float, None]]:
    """
    Creates a placeholder metrics dictionary for cases where a network
    is unavailable or too small.

    Returns
    -------
    Dict[str, Union[str, int, float, None]]
        A dictionary containing placeholder metrics with a message
        for unavailable or small networks.
    """
    return {
        "Threshold": threshold,
        "NetworkType": network_type,
        "Dataset": dataset_name,
        "Nodes": n_nodes if n_nodes is not None else message,
        "Edges": n_edges if n_edges is not None else message,
        "Number_communities": message,
        "Modularity": message,
        "Avg_degree": message,
        "Avg_clustering_coefficient": message,
        "Density": message,
    }


def process_and_analyze_networks(
    base_dir: str, network_type: str, output_dir: str, weight_column: str
) -> Tuple[pd.DataFrame, Union[nx.Graph, None], Union[str, None]]:
    """
    Processes all network folders in a base directory, calculates metrics,
    saves them to a CSV file, and identifies the best network based on
    modularity and clustering coefficient.

    Parameters
    ----------
    base_dir : str
        The base directory containing network subfolders (e.g., 'cclasso_0.10').
    network_type : str
        The type of network to process (e.g., 'cclasso', 'spring'). Used for folder matching.
    output_dir : str
        The directory where the analysis results (metrics CSV and best network clustering) will be saved.
    weight_column : str
        The column to use as edge weight for network metrics ('adja', 'asso', or 'diss').

    Returns
    -------
    Tuple[pd.DataFrame, Union[nx.Graph, None], Union[str, None]]
        A tuple containing:
        - metrics_df (pd.DataFrame): A DataFrame with all calculated network metrics.
        - best_network_graph (nx.Graph or None): The NetworkX graph object identified as the best, or
          None if no best network is found.
        - best_threshold (str or None): The threshold of the best network, or None.
    """
    metrics_list = []
    top_networks_by_modularity = []

    # The "dataset" is the base directory itself
    dataset_name = os.path.basename(os.path.normpath(base_dir))

    # Process each network folder within the base directory
    for folder_name in os.listdir(base_dir):
        if not folder_name.startswith(network_type):
            continue

        network_folder_path = os.path.join(base_dir, folder_name)
        if not os.path.isdir(network_folder_path):
            continue

        try:
            threshold: str = folder_name.split("_")[-1]
            threshold_numeric = float(threshold)
        except ValueError:
            print(
                f"Warning: Could not extract numeric threshold from folder name: {folder_name}. Skipping."
            )
            continue

        network_files = [
            f
            for f in os.listdir(network_folder_path)
            if f.endswith(".csv") and "nosinglt" in f
        ]

        if not network_files:
            # Case where no network file is available
            metrics_list.append(
                create_metrics_inconsistent(
                    threshold, network_type, dataset_name, "No network available"
                )
            )
            continue

        # Keep only the first file
        file = network_files[0]
        file_path = os.path.join(network_folder_path, file)
        print(f"Analyzing network: {file_path}")

        # Create a graph from the file
        G = None
        try:
            G = create_networkx_graph(file_path, weight_column)
        except ValueError as e:
            print(f"Error creating graph from {file_path}: {e}. Skipping.")
            metrics_list.append(
                create_metrics_inconsistent(
                    threshold, network_type, dataset_name, f"Error: {e}"
                )
            )
            continue
        except Exception as e:
            print(
                f"An unexpected error occurred creating graph from {file_path}: {e}. Skipping."
            )
            metrics_list.append(
                create_metrics_inconsistent(
                    threshold, network_type, dataset_name, f"Unexpected Error: {e}"
                )
            )
            continue

        # Case where network is too small
        if G is None or G.number_of_edges() <= 3:
            metrics_list.append(
                create_metrics_inconsistent(
                    threshold,
                    network_type,
                    dataset_name,
                    "Network too small (<= 3 edges)",
                    G.number_of_nodes() if G else 0,
                    G.number_of_edges() if G else 0,
                )
            )
            continue

        # Ensure weights are absolute for metrics calculation, as some metrics do not support signed values.
        # This transformation applies to 'asso' (converting negative to positive) and has no effect on
        # 'adja' or 'diss' if they are already non-negative.
        print(
            f"Transforming '{weight_column}' weights to absolute values for metrics calculation."
        )
        G_for_metrics = transform_weights_to_absolute(
            G.copy(), weight_attribute=weight_column
        )

        # Calculate network metrics and append to the list
        network_metrics = calculate_network_metrics(G_for_metrics, weight_column)
        metrics = {
            "Threshold": threshold,
            "NetworkType": network_type,
            "Dataset": dataset_name,
            **network_metrics,
        }
        metrics_list.append(metrics)
        print(f"Network metrics calculated for {file_path}")

        # Track top networks only if modularity is a positive number
        modularity_value = network_metrics["Modularity"]
        if isinstance(modularity_value, (int, float)) and modularity_value > 0:
            top_networks_by_modularity.append(
                (
                    float(modularity_value),  # Ensure it's float for sorting
                    float(
                        network_metrics["Avg_clustering_coefficient"]
                    ),  # Ensure it's float for sorting
                    G,
                    threshold,
                )
            )
            top_networks_by_modularity.sort(key=lambda x: x[0], reverse=True)
            top_networks_by_modularity = top_networks_by_modularity[:3]

    # Create and save the DataFrame with all metrics
    metrics_df = pd.DataFrame(metrics_list)
    output_csv_path = os.path.join(
        output_dir, f"{dataset_name}_{network_type}_{weight_column}_metrics.csv"
    )
    metrics_df.to_csv(output_csv_path, index=False)
    print(f"Unified network metrics saved to {output_csv_path}")

    # Select the best network from the top three based on the highest clustering coefficient
    best_network_graph, best_threshold = None, None
    if top_networks_by_modularity:
        # max() on a list of tuples compares element by element; here we sort by the second element (clustering coeff)
        best_network = max(top_networks_by_modularity, key=lambda x: x[1])
        _, _, best_network_graph, best_threshold = best_network

    return metrics_df, best_network_graph, best_threshold


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Perform topological analysis on microbial association networks genarated with the NetComi package with different thresholds."
    )
    parser.add_argument(
        "--base_dir",
        type=str,
        help="Base directory containing network folders (e.g., 'cclasso_0.10', 'cclasso_0.15').",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        help="Directory to save the analysis results (metrics CSV and best network clustering).",
    )
    parser.add_argument(
        "--network_type",
        type=str,
        default="cclasso",
        help="Type of network to process (e.g., 'cclasso', 'spring'). Used for folder matching.",
    )
    parser.add_argument(
        "--weight_column",
        type=str,
        default="adja",
        choices=["adja", "asso", "diss"],
        help="Column to use as edge weight for network metrics (e.g., 'adja', 'asso', 'diss').",
    )

    args = parser.parse_args()

    # Process the networks, get metrics and the best network
    metrics_df, best_network_graph, best_threshold = process_and_analyze_networks(
        args.base_dir, args.network_type, args.output_dir, args.weight_column
    )

    # If a best network was found, export its node clustering information
    if best_network_graph:
        dataset_name: str = os.path.basename(os.path.normpath(args.base_dir))
        print(
            f"\nBest network found at threshold: {best_threshold} for dataset: {dataset_name}"
        )

        # Ensure weights are absolute for best network clustering, as some metrics do not support signed values.
        # This transformation applies to 'asso' (converting negative to positive) and has no effect on
        # 'adja' or 'diss' if they are already non-negative.
        print(
            f"    Transforming '{args.weight_column}' weights to absolute values for best network clustering."
        )
        best_network_graph_for_clustering = transform_weights_to_absolute(
            best_network_graph.copy(), weight_attribute=args.weight_column
        )

        # Calculate communities for the best network
        communities: List[List[Any]] = nx.algorithms.community.louvain_communities(
            best_network_graph_for_clustering, weight=args.weight_column, seed=12
        )
        node_community_dict: Dict[Any, int] = {
            node: cid for cid, community in enumerate(communities) for node in community
        }
        df_node_community: pd.DataFrame = pd.DataFrame(
            node_community_dict.items(), columns=["Node", "Community"]
        )

        # Define and create the output directory for clustering results
        best_net_output_dir: str = os.path.join(args.output_dir, "Best_nets/Clustering")
        os.makedirs(best_net_output_dir, exist_ok=True)

        # Save the node clustering data to a CSV file
        clustering_output_path: str = os.path.join(
            best_net_output_dir,
            f"{dataset_name}_best_net_{best_threshold}_node_clustering_{args.weight_column}.csv",
        )
        df_node_community.to_csv(clustering_output_path, index=False)
        print(
            f"Node clustering information for the best network saved to {clustering_output_path}"
        )
    else:
        print("\nNo single best network could be determined based on the criteria.")
