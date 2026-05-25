"""
preprocessing.py
HVG-selectie → scaling → PCA → BBKNN batch-correctie.
"""

import scanpy as sc
import anndata as ad


def run_preprocessing(adata: ad.AnnData, n_top_genes: int = 2000, n_comps: int = 20) -> ad.AnnData:
    """
    Voer de volledige preprocessing pipeline uit:
    HVG-selectie, scaling, PCA en BBKNN batch-correctie.

    Parameters
    ----------
    adata       : gecombineerd AnnData object (output van data_loader)
    n_top_genes : aantal highly variable genes om te selecteren (default 2000)
    n_comps     : aantal PCA componenten (default 20)

    Returns
    -------
    adata : preprocessed AnnData object
    """
    # --- 1. HVG-selectie ---
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        flavor="seurat",
        batch_key="sample"
    )

    # --- Raw opslaan VOOR HVG-filter (nodig voor CellTypist) ---
    adata.raw = adata

    # --- 2. Filter naar alleen HVGs ---
    adata = adata[:, adata.var.highly_variable].copy()

    # --- 3. Schalen & PCA ---
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=n_comps)

    # --- 4. BBKNN batch-correctie ---
    sc.external.pp.bbknn(
        adata,
        batch_key="sample",
        neighbors_within_batch=9,
        n_pcs=4
    )

    print("Preprocessing voltooid.")
    return adata
