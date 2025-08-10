# %% Imports
import pandas as pd
import networkx as nx
import numpy as np
from scipy.stats import lognorm
import matplotlib.pyplot as plt

# %% Load edge lists
# Parameters
lnormFit = True  # Set to False to use empirical threshold

# Load edge list
df_net_gen = pd.read_csv(
    "../../Output/Best_nets/PGE_and_NOPGE_combined_gen/cclasso_nosinglt_70_edgelist.csv",
    sep=",",
)
df_net_sp = pd.read_csv(
    "../../Output/Best_nets/PGE_and_NOPGE_combined_sp/cclasso_nosinglt_65_edgelist.csv",
    sep=",",
)

# %% Build the networks using 'adja' as weight
net_gen = nx.Graph()
for _, row in df_net_gen.iterrows():
    net_gen.add_edge(row["v1"], row["v2"], weight=row["adja"])

net_sp = nx.Graph()
for _, row in df_net_sp.iterrows():
    net_sp.add_edge(row["v1"], row["v2"], weight=row["adja"])


# %% Compute centrality metrics
def calc_centrality_measures(G):
    centralities = {}

    # Degree and betweenness work on disconnected graphs
    centralities["degree"] = nx.degree_centrality(G)
    centralities["betweenness"] = nx.betweenness_centrality(G, weight="weight")

    if not nx.is_connected(G):
        print(
            "Graph is not connected. Using largest connected component for closeness and eigenvector."
        )
        largest_cc = max(nx.connected_components(G), key=len)
        G = G.subgraph(largest_cc).copy()

    # These require connected components
    centralities["closeness"] = nx.closeness_centrality(G)
    centralities["eigenvector"] = nx.eigenvector_centrality(
        G, weight="weight", max_iter=1000
    )

    return centralities


centralities_net_gen = calc_centrality_measures(net_gen)
centralities_net_sp = calc_centrality_measures(net_sp)
# %% Identify hubs
hubs_net_gen = {}

for name, centrality in centralities_net_gen.items():
    values = np.array(list(centrality.values()))

    if lnormFit:
        # Filter out zero or negative values
        values_pos = values[values > 0]
        if len(values_pos) < 2:
            print(
                f"Not enough positive values to fit log-normal for {name}. Using empirical threshold instead."
            )
            threshold = np.percentile(values, 95)
        else:
            shape, loc, scale = lognorm.fit(values_pos, floc=0)
            threshold = lognorm.ppf(0.95, shape, loc=loc, scale=scale)
    else:
        # Empirical 95th percentile
        threshold = np.percentile(values, 95)

    # Nodes above threshold
    hubs_net_gen[name] = [node for node, val in centrality.items() if val > threshold]
    print(f"\n>>> Hubs by {name} centrality (threshold = {threshold:.4f}):")
    for node in hubs_net_gen[name]:
        print(f" - {node}: {centrality[node]:.4f}")

# %%

hubs_net_sp = {}

for name, centrality in centralities_net_sp.items():
    values = np.array(list(centrality.values()))

    if lnormFit:
        # Filter out zero or negative values
        values_pos = values[values > 0]
        if len(values_pos) < 2:
            print(
                f"Not enough positive values to fit log-normal for {name}. Using empirical threshold instead."
            )
            threshold = np.percentile(values, 95)
        else:
            shape, loc, scale = lognorm.fit(values_pos, floc=0)
            threshold = lognorm.ppf(0.85, shape, loc=loc, scale=scale)
    else:
        # Empirical 95th percentile
        threshold = np.percentile(values, 95)

    # Nodes above threshold
    hubs_net_sp[name] = [node for node, val in centrality.items() if val > threshold]
    print(f"\n>>> Hubs by {name} centrality (threshold = {threshold:.4f}):")
    for node in hubs_net_sp[name]:
        print(f" - {node}: {centrality[node]:.4f}")

# %%
