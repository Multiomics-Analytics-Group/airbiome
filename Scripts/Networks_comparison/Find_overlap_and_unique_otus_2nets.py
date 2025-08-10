# %%
import pandas as pd


def extract_otus(network_path):
    """
    Extract OTUs from a CSV network file.
    """
    df = pd.read_csv(network_path)
    return set(df["v1"]).union(df["v2"])


def compare_networks(
    net_file1,
    net_file2,
    label1="Network1",
    label2="Network2",
    overlap_csv_path=None,
    unique_net1_csv_path=None,
    unique_net2_csv_path=None,
):
    """
    Compare two networks by OTUs and report overlaps.
    Optionally save the overlapping OTUs to a CSV file.
    """
    otus1 = extract_otus(net_file1)
    otus2 = extract_otus(net_file2)
    common_otus = otus1 & otus2
    unique_otus1 = otus1 - otus2
    unique_otus2 = otus2 - otus1

    print(f"{label1} OTUs: {len(otus1)}")
    print(f"{label2} OTUs: {len(otus2)}")
    print(f"Overlapping OTUs: {len(common_otus)}")
    print(f"{label1} unique OTUs: {len(unique_otus1)}")
    print(f"{label2} unique OTUs: {len(unique_otus2)}")

    # Optional metrics
    union_size = len(otus1 | otus2)
    jaccard = len(common_otus) / union_size if union_size > 0 else 0
    print(f"Jaccard Similarity: {jaccard:.4f}")
    print(f"Percent overlap (vs {label1}): {100 * len(common_otus) / len(otus1):.2f}%")
    print(f"Percent overlap (vs {label2}): {100 * len(common_otus) / len(otus2):.2f}%")

    # Save overlapping OTUs
    if overlap_csv_path:
        pd.DataFrame({"overlapping_otus": list(common_otus)}).to_csv(
            overlap_csv_path, index=False
        )
        print(f"Overlapping OTUs saved to: {overlap_csv_path}")

    # Save unique OTUs for network1
    if unique_net1_csv_path:
        pd.DataFrame({"unique_otus": sorted(unique_otus1)}).to_csv(
            unique_net1_csv_path, index=False
        )
        print(f"Unique {label1} OTUs saved to: {unique_net1_csv_path}")

    # Save unique OTUs for network2
    if unique_net2_csv_path:
        pd.DataFrame({"unique_otus": sorted(unique_otus2)}).to_csv(
            unique_net2_csv_path, index=False
        )
        print(f"Unique {label2} OTUs saved to: {unique_net2_csv_path}")


# %%
PGE_net = "../../Output/Network_inference/PGE_and_NOPGE_05/PGE/cclasso_0.50/PGE_nosinglt_cclasso_50_edgelist.csv"
NOPGE_net = "../../Output/Network_inference/PGE_and_NOPGE_05/NOPGE/cclasso_0.50/NOPGE_nosinglt_cclasso_50_edgelist.csv"
overlap_otus_path = (
    "../../Output/Best_nets/otu_overlap_airbiome_cclasso_05_nosinglet.csv"
)
unique_otus_pge_path = "../../Output/Best_nets/unique_pge_otus_cclasso_05_nosinglet.csv"

compare_networks(
    PGE_net,
    NOPGE_net,
    label1="PGE",
    label2="NO_PGE",
    overlap_csv_path=overlap_otus_path,
    unique_net1_csv_path=unique_otus_pge_path,
)

# %%
