"""
training.py
Hoofdscript dat de volledige pipeline uitvoert:
data laden → preprocessing → clustering → annotatie → model training.

Gebruik:
    python training.py
"""

import scanpy as sc
from data_loader   import load_and_merge
from preprocessing import run_preprocessing
from dimensionality import run_clustering, plot_umap, plot_cluster_composition
from annotation    import compute_marker_genes, annotate_clusters, train_celltypist_model

# --- Instellingen ---
sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=100, frameon=False)

BF_PATH    = "../data/BF_count_mat.h5ad"
P_PATH     = "../data/P_count_mat.h5ad"
MODEL_PATH = "Staph_Aureus_BF_P_Model.pkl"


def main():
    # 1. Data laden en samenvoegen
    adata = load_and_merge(BF_PATH, P_PATH)

    # 2. Preprocessing: HVG → scaling → PCA → BBKNN
    adata = run_preprocessing(adata)

    # 3. Dimensionaliteitsreductie en clustering
    adata = run_clustering(adata)
    plot_umap(adata)
    dist_norm = plot_cluster_composition(adata)
    print("\nCluster-compositie (% Biofilm vs Planktonic):")
    print(dist_norm.round(2))

    # 4. Markergen-analyse
    marker_table = compute_marker_genes(adata)

    # 5. Cluster-annotatie
    # LET OP: bekijk de markergen-tabel en de compositie-grafiek
    # en pas DEFAULT_CLUSTER_NAMES in annotation.py aan indien nodig.
    adata = annotate_clusters(adata)

    # 6. CellTypist model trainen en opslaan
    train_celltypist_model(adata, model_path=MODEL_PATH)


if __name__ == "__main__":
    main()
