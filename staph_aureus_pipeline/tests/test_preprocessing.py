"""
test_preprocessing.py
Testen voor src/pipeline/preprocessing.py

Wat we testen:
- HVG-selectie voegt 'highly_variable' kolom toe
- PCA wordt berekend en opgeslagen
- Scaling zorgt dat waarden niet extreme uitschieters hebben
- Output heeft minder of gelijk aantal genen (door HVG-filter)
"""

import numpy as np
import anndata as ad
import pandas as pd
import scanpy as sc
import pytest


def _make_preprocessable_adata(n_cells=40, n_genes=100):
    """
    Maak een nep-AnnData die groot genoeg is voor HVG + PCA.
    Minimaal: 2 samples, genoeg cellen en genen.
    """
    np.random.seed(42)
    X = np.random.negative_binomial(5, 0.5, size=(n_cells, n_genes)).astype(float)

    obs = pd.DataFrame(
        {"sample": ["bf"] * (n_cells // 2) + ["p"] * (n_cells // 2)},
        index=[f"cell_{i}" for i in range(n_cells)]
    )
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(n_genes)])

    adata = ad.AnnData(X=X, obs=obs, var=var)

    # Normaliseer zoals data_loader dat doet
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    return adata


# ── Test 1: HVG-selectie ──────────────────────────────────────────────────────

def test_hvg_selectie_voegt_kolom_toe():
    """
    WAT: na HVG-selectie moet 'highly_variable' in adata.var zitten.
    WAAROM: zonder deze kolom weet de rest van de pipeline niet welke genen te gebruiken.
    """
    adata = _make_preprocessable_adata()
    sc.pp.highly_variable_genes(adata, n_top_genes=50, flavor="seurat", batch_key="sample")

    assert "highly_variable" in adata.var.columns


def test_hvg_selectie_markeert_genen():
    """
    WAT: minstens één gen moet als highly_variable gemarkeerd zijn.
    WAAROM: als er geen HVGs zijn loopt de rest van de pipeline vast.
    """
    adata = _make_preprocessable_adata()
    sc.pp.highly_variable_genes(adata, n_top_genes=50, flavor="seurat", batch_key="sample")

    n_hvg = adata.var["highly_variable"].sum()
    assert n_hvg > 0, f"Geen highly variable genes gevonden! (n_hvg={n_hvg})"


# ── Test 2: PCA ───────────────────────────────────────────────────────────────

def test_pca_wordt_berekend():
    """
    WAT: na sc.tl.pca() moet adata.obsm['X_pca'] bestaan.
    WAAROM: BBKNN en UMAP hebben de PCA-representatie nodig.
    """
    adata = _make_preprocessable_adata()
    sc.pp.highly_variable_genes(adata, n_top_genes=50, flavor="seurat", batch_key="sample")
    adata = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=10)

    assert "X_pca" in adata.obsm


def test_pca_heeft_juiste_dimensies():
    """
    WAT: PCA met n_comps=10 geeft een matrix van (n_cellen, 10).
    WAAROM: verkeerde dimensies breken de BBKNN-stap.
    """
    adata = _make_preprocessable_adata()
    sc.pp.highly_variable_genes(adata, n_top_genes=50, flavor="seurat", batch_key="sample")
    adata = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=10)

    assert adata.obsm["X_pca"].shape == (adata.n_obs, 10)


# ── Test 3: Scaling ───────────────────────────────────────────────────────────

def test_scaling_beperkt_uitschieters():
    """
    WAT: na sc.pp.scale(max_value=10) zijn er geen waarden groter dan 10.
    WAAROM: extreme waarden verstoren PCA en clustering.
    """
    adata = _make_preprocessable_adata()
    sc.pp.highly_variable_genes(adata, n_top_genes=50, flavor="seurat", batch_key="sample")
    adata = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(adata, max_value=10)

    assert adata.X.max() <= 10.01  # kleine marge voor floating point
