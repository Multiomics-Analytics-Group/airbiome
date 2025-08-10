# %%
import pandas as pd


def load_network(file_path):
    """
    Load a network file and return a set of undirected edges and original DataFrame.
    """
    df = pd.read_csv(file_path)[["v1", "v2"]]
    edges = set(map(frozenset, df.values))  # Undirected edges
    return edges, df


def compare_two_networks(edges1, edges2, label1, label2):
    """
    Compare two edge sets and return overlap stats and common edges.
    """
    intersection = edges1 & edges2
    union = edges1 | edges2

    jaccard_similarity = len(intersection) / len(union) if union else 0
    percent_overlap_1 = len(intersection) / len(edges1) * 100 if edges1 else 0
    percent_overlap_2 = len(intersection) / len(edges2) * 100 if edges2 else 0

    print(f"Comparison of networks: {label1} vs {label2}")
    print(f"{label1} edges: {len(edges1)}")
    print(f"{label2} edges: {len(edges2)}")
    print(f"Overlapping edges: {len(intersection)}")
    print(f"Jaccard similarity: {jaccard_similarity:.4f}")
    print(f"Percent overlap (vs {label1}): {percent_overlap_1:.2f}%")
    print(f"Percent overlap (vs {label2}): {percent_overlap_2:.2f}%")

    return intersection


def save_overlapping_edges(intersection_set, output_path):
    """
    Save the overlapping edges to a CSV file.
    """
    # Convert frozensets to sorted edge pairs (v1, v2)
    edge_list = [sorted(list(edge)) for edge in intersection_set]
    df = pd.DataFrame(edge_list, columns=["v1", "v2"])
    df.to_csv(output_path, index=False)
    print(f"Overlapping edges saved to {output_path}")


# %% Main execution
PGE_net = "../../Output/Network_inference/PGE_and_NOPGE_05/PGE/cclasso_0.50/PGE_nosinglt_cclasso_50_edgelist.csv"
NOPGE_net = "../../Output/Network_inference/PGE_and_NOPGE_05/NOPGE/cclasso_0.50/NOPGE_nosinglt_cclasso_50_edgelist.csv"
overlap_edges_path = "../../Output/edge_overlap_airbiome_cclasso_05_nosinglet.csv"

PGE_edges, _ = load_network(PGE_net)
NOPGE_edges, _ = load_network(NOPGE_net)

overlapping_edges = compare_two_networks(PGE_edges, NOPGE_edges, "PGE", "NO_PGE")
save_overlapping_edges(overlapping_edges, overlap_edges_path)
# %%
