"""
annotation.py
Markergen-analyse, handmatige cluster-annotatie en CellTypist model training.
"""

import scanpy as sc
import anndata as ad
import pandas as pd
import celltypist


# Standaard cluster-namen op basis van het Korshoj & Kielian (2024) paper.
# Pas deze aan op basis van jouw markergen-tabel en cluster-compositie.
DEFAULT_CLUSTER_NAMES = {
    "0": "Transitional",                     # ~50/50 BF/P
    "1": "Biofilm_Transcriptionally_Active",  # citB, ltaS, isdH
    "2": "Planktonic_Active",                 # rpoB, rpsC
    "3": "Biofilm_Virulence",                 # clfB, fnbB, arcA
    "4": "Planktonic_Stationary",             # qoxA, hemY
    "5": "Biofilm_Stress_Metabolism",         # gpmA, polX, fdaB
    "6": "Biofilm_Replication"                # nrdE
}


def compute_marker_genes(adata: ad.AnnData) -> pd.DataFrame:
    """
    Bereken markergen per cluster (Wilcoxon) en print top 5.

    Returns
    -------
    marker_table : DataFrame met top markergenen per cluster
    """
    sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
    sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False)

    marker_table = pd.DataFrame(
        adata.uns["rank_genes_groups"]["names"]
    ).head(5)

    print("Top markergen per cluster:")
    print(marker_table)
    return marker_table


def annotate_clusters(
    adata: ad.AnnData,
    cluster_names: dict = None
) -> ad.AnnData:
    """
    Wijs celtype-namen toe aan Leiden-clusters.

    Parameters
    ----------
    adata         : AnnData met 'leiden' kolom
    cluster_names : dict van cluster-ID (str) naar celtypnaam.
                    Standaard: DEFAULT_CLUSTER_NAMES uit dit bestand.

    Returns
    -------
    adata : AnnData met extra 'cell_type' kolom
    """
    if cluster_names is None:
        cluster_names = DEFAULT_CLUSTER_NAMES

    # Controleer of alle clusters gedekt zijn
    clusters_in_data = sorted(adata.obs["leiden"].unique())
    clusters_in_dict = sorted(cluster_names.keys())
    print(f"Clusters in data:  {clusters_in_data}")
    print(f"Clusters in dict:  {clusters_in_dict}")

    adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_names)

    n_unmapped = adata.obs["cell_type"].isna().sum()
    print(f"Niet-gemapte cellen: {n_unmapped}")
    if n_unmapped > 0:
        print("Waarschuwing: pas DEFAULT_CLUSTER_NAMES aan voor de ontbrekende clusters.")

    sc.pl.umap(
        adata,
        color="cell_type",
        legend_loc="on data",
        title="Cell type annotation"
    )
    return adata


def train_celltypist_model(
    adata: ad.AnnData,
    model_path: str = "Staph_Aureus_BF_P_Model.pkl",
    top_genes: int = 300
) -> celltypist.Model:
    """
    Train een CellTypist model op basis van de annotaties in 'cell_type'.

    Parameters
    ----------
    adata      : AnnData met 'cell_type' kolom (geen NaN)
    model_path : pad om het model op te slaan
    top_genes  : aantal top genen voor feature-selectie

    Returns
    -------
    model : getraind CellTypist model
    """
    # Veiligheidscheck
    n_missing = adata.obs["cell_type"].isna().sum()
    assert n_missing == 0, f"Los {n_missing} ontbrekende labels op voor training!"

    model = celltypist.train(
        adata,
        labels="cell_type",
        feature_selection=True,
        top_genes=top_genes
    )

    model.write(model_path)
    print(f"Model opgeslagen als: {model_path}")
    print("\nTop features per celtype:")
    print(model.features)

    return model
