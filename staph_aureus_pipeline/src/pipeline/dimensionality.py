"""
dimensionality.py
UMAP-embedding, Leiden-clustering en visualisaties (UMAP + cluster compositie).
"""

import scanpy as sc
import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt


def run_clustering(
    adata: ad.AnnData,
    resolution: float = 0.205,
    min_dist: float = 0.24,
    spread: float = 0.21
) -> ad.AnnData:
    """
    Bereken UMAP en Leiden-clustering.

    Parameters
    ----------
    adata      : preprocessed AnnData (output van preprocessing)
    resolution : Leiden resolutie (0.205 geeft ~7 clusters zoals in het paper)
    min_dist   : UMAP min_dist parameter
    spread     : UMAP spread parameter

    Returns
    -------
    adata : AnnData met 'leiden' en 'Cell identity' kolommen
    """
    # --- UMAP & Leiden ---
    sc.tl.umap(adata, min_dist=min_dist, spread=spread)
    sc.tl.leiden(adata, resolution=resolution)

    n_clusters = adata.obs["leiden"].nunique()
    print(f"Aantal gevonden clusters: {n_clusters} (verwacht: 7)")

    # --- Conditie-labels toevoegen ---
    actual_labels = adata.obs["sample"].unique()
    label_map = {actual_labels[0]: "Biofilm", actual_labels[1]: "Planktonic"}
    adata.obs["Cell identity"] = adata.obs["sample"].map(label_map)

    return adata


def plot_umap(adata: ad.AnnData) -> None:
    """Plot UMAP gekleurd op conditie en Leiden-cluster."""
    sc.pl.umap(
        adata,
        color=["Cell identity", "leiden"],
        title=["Condition", "Leiden clusters"],
        wspace=0.4,
        frameon=False
    )


def plot_cluster_composition(adata: ad.AnnData) -> pd.DataFrame:
    """
    Maak een gestapeld staafdiagram van cluster-compositie (Figure 2D).

    Returns
    -------
    dist_norm : genormaliseerde proportie-tabel (clusters × celtype)
    """
    dist = pd.crosstab(adata.obs["leiden"], adata.obs["Cell identity"])
    dist_norm = dist.div(dist.sum(axis=1), axis=0)

    color_map = {"Biofilm": "steelblue", "Planktonic": "tomato"}
    fig, axes = plt.subplots(1, 2, figsize=(14, 4))

    # Absoluut
    dist.plot(kind="bar", ax=axes[0], color=color_map, stacked=True)
    axes[0].set_ylabel("Number of cells")
    axes[0].set_xlabel("Leiden Cluster")
    axes[0].set_title("Absolute cell types per cluster")
    axes[0].tick_params(axis="x", rotation=0)
    axes[0].legend(title="Condition", frameon=False)

    # Genormaliseerd
    dist_norm.plot(kind="bar", ax=axes[1], color=color_map, stacked=True)
    axes[1].set_ylabel("Proportion of cells")
    axes[1].set_xlabel("Leiden Cluster")
    axes[1].set_title("Proportional cell types per cluster")
    axes[1].tick_params(axis="x", rotation=0)
    axes[1].axhline(0.5, color="black", linestyle="--", linewidth=0.8)
    axes[1].legend(
        title="Condition", frameon=False,
        loc="center left", bbox_to_anchor=(1, 0.5)
    )

    plt.tight_layout()
    plt.show()

    return dist_norm
