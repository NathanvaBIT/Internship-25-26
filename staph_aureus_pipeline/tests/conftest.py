"""
conftest.py
Gedeelde nep-data (fixtures) die alle testbestanden kunnen gebruiken.
Pytest laadt dit bestand automatisch — je hoeft het nooit zelf te importeren.
"""

import pytest
import numpy as np
import anndata as ad
import pandas as pd
import scipy.sparse as sp


def _make_adata(n_cells: int, n_genes: int, label: str) -> ad.AnnData:
    """Hulpfunctie: maak een klein nep-AnnData object."""
    # Willekeurige count-matrix (gehele getallen, zoals echte scRNA-seq data)
    X = np.random.randint(0, 50, size=(n_cells, n_genes)).astype(float)

    # Gennamen: een paar ribosomale genen ertussen zodat de filter testbaar is
    gene_names = [f"gene_{i}" for i in range(n_genes)]
    gene_names[0] = "rpsA"   # ribosomaal → moet gefilterd worden
    gene_names[1] = "rplB"   # ribosomaal → moet gefilterd worden

    obs = pd.DataFrame(index=[f"{label}_cell_{i}" for i in range(n_cells)])
    var = pd.DataFrame(index=gene_names)

    return ad.AnnData(X=X, obs=obs, var=var)


@pytest.fixture
def adata_bf():
    """Klein nep BF (Biofilm) dataset — 20 cellen, 50 genen."""
    adata = _make_adata(n_cells=20, n_genes=50, label="bf")
    # Zorg dat counts hoog genoeg zijn om de count-filter (>= 7) te passeren
    adata.X[:, 2:] = np.random.randint(10, 50, size=(20, 48))
    return adata


@pytest.fixture
def adata_p():
    """Klein nep P (Planktonic) dataset — 20 cellen, 50 genen."""
    adata = _make_adata(n_cells=20, n_genes=50, label="p")
    # Zorg dat counts hoog genoeg zijn om de count-filter (>= 28) te passeren
    adata.X[:, 2:] = np.random.randint(30, 50, size=(20, 48))
    return adata


@pytest.fixture
def adata_merged(adata_bf, adata_p):
    """Gecombineerd nep-dataset zoals load_and_merge() dat teruggeeft."""
    import anndata as ad
    import scanpy as sc

    # Ribo-filter nabootsen
    for adata in [adata_bf, adata_p]:
        adata.var_names_make_unique()
        adata.var["ribo"] = adata.var_names.str.lower().str.startswith(
            ("rps", "rpl", "rpm", "rrs", "rrl", "rrf")
        )
        adata = adata[:, ~adata.var["ribo"]].copy()

    # Normaliseer
    sc.pp.normalize_total(adata_bf, target_sum=1e4)
    sc.pp.log1p(adata_bf)
    sc.pp.normalize_total(adata_p, target_sum=1e4)
    sc.pp.log1p(adata_p)

    merged = ad.concat(
        {"bf": adata_bf, "p": adata_p},
        label="sample",
        join="inner"
    )
    merged.obs_names_make_unique()
    return merged
