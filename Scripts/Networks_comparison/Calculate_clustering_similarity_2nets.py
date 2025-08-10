# %%
import pandas as pd
import networkx as nx
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics import normalized_mutual_info_score


def load_network_from_csv(filepath):
    """
    Load an edge list CSV and return a NetworkX graph.
    """
    df = pd.read_csv(filepath)
    G = nx.Graph()
    G.add_edges_from(df[["v1", "v2"]].values)
    return G


def compute_louvain_communities(G):
    """
    Compute Louvain communities using NetworkX.
    Returns a dict: node -> community_id
    """
    communities = nx.algorithms.community.louvain_communities(G, seed=12)
    return {
        node: cid for cid, community in enumerate(communities) for node in community
    }


def compute_ari(comm1, comm2):
    shared_nodes = list(set(comm1.keys()).intersection(set(comm2.keys())))
    labels1 = [comm1[n] for n in shared_nodes]
    labels2 = [comm2[n] for n in shared_nodes]
    ari = adjusted_rand_score(labels1, labels2)
    return ari


def compute_nmi(comm1, comm2):
    shared_nodes = list(set(comm1.keys()).intersection(set(comm2.keys())))
    labels1 = [comm1[n] for n in shared_nodes]
    labels2 = [comm2[n] for n in shared_nodes]
    nmi = normalized_mutual_info_score(labels1, labels2)
    return nmi


# %% Main execution
PGE_net = "../../Output/Network_inference/PGE_and_NOPGE_05/PGE/cclasso_0.50/PGE_net_cclasso_50_edgelist.csv"
NOPGE_net = "../../Output/Network_inference/PGE_and_NOPGE_05/NOPGE/cclasso_0.50/NOPGE_net_cclasso_50_edgelist.csv"
overlap_comm_path = "../../Output/comm_overlap_airbiome_cclasso_05.csv"

# Load and compute communities
G1 = load_network_from_csv(PGE_net)
G2 = load_network_from_csv(NOPGE_net)

comm1 = compute_louvain_communities(G1)
comm2 = compute_louvain_communities(G2)

# Compare and save results
ari = compute_ari(comm1, comm2)
nmi = compute_nmi(comm1, comm2)

# Print the scores
print(f"Adjusted Rand Index (ARI): {ari:.4f}")
print(f"Normalized Mutual Information (NMI): {nmi:.4f}")

# %%
