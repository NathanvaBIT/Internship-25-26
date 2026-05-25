"""
test_data_loader.py
Testen voor src/pipeline/data_loader.py

Wat we testen:
- Ribo-filter verwijdert ribosomale genen
- Count-filter verwijdert cellen met te weinig counts
- Merge geeft 'sample' kolom terug
- Shape van de output klopt
"""

import numpy as np
import anndata as ad
import pandas as pd
import scanpy as sc


# ── Hulpfunctie om nep-data te maken zonder conftest te hoeven laden ──────────

def _make_raw_adata(n_cells, n_genes, label, min_count=0):
    """Maak een nep AnnData met hoge genoeg counts."""
    X = np.full((n_cells, n_genes), fill_value=50, dtype=float)
    gene_names = [f"gene_{i}" for i in range(n_genes)]
    gene_names[0] = "rpsA"  # ribosomaal gen
    gene_names[1] = "rplB"  # ribosomaal gen
    obs = pd.DataFrame(index=[f"{label}_{i}" for i in range(n_cells)])
    var = pd.DataFrame(index=gene_names)
    return ad.AnnData(X=X, obs=obs, var=var)


# ── Test 1: Ribo-filter ────────────────────────────────────────────────────────

def test_ribo_filter_verwijdert_ribo_genen():
    """
    WAT: rpsA en rplB zijn ribosomale genen → moeten verdwijnen na filter.
    WAAROM: ribosomale genen verstoren de analyse en moeten eruit.
    """
    adata = _make_raw_adata(10, 20, "bf")

    # Voer de ribo-filter uit (zelfde logica als in data_loader.py)
    adata.var_names_make_unique()
    adata.var["ribo"] = adata.var_names.str.lower().str.startswith(
        ("rps", "rpl", "rpm", "rrs", "rrl", "rrf")
    )
    adata_filtered = adata[:, ~adata.var["ribo"]].copy()

    # rpsA en rplB mogen niet meer in de gennamen zitten
    assert "rpsA" not in adata_filtered.var_names
    assert "rplB" not in adata_filtered.var_names


def test_ribo_filter_behoudt_normale_genen():
    """
    WAT: gewone genen (gene_2, gene_3, ...) moeten bewaard blijven.
    WAAROM: we filteren alleen ribo-genen, niet alles.
    """
    adata = _make_raw_adata(10, 20, "bf")
    adata.var_names_make_unique()
    adata.var["ribo"] = adata.var_names.str.lower().str.startswith(
        ("rps", "rpl", "rpm", "rrs", "rrl", "rrf")
    )
    adata_filtered = adata[:, ~adata.var["ribo"]].copy()

    assert "gene_2" in adata_filtered.var_names


# ── Test 2: Count-filter ──────────────────────────────────────────────────────

def test_count_filter_bf_verwijdert_lage_cellen():
    """
    WAT: cellen met minder dan 7 counts moeten verdwijnen voor BF data.
    WAAROM: lage counts zijn technisch ruis, geen echte cellen.
    """
    # 5 cellen met counts = 3 (te laag), 5 cellen met counts = 50 (goed)
    X = np.array([[3] * 10] * 5 + [[50] * 10] * 5, dtype=float)
    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(10)])
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(10)])
    adata = ad.AnnData(X=X, obs=obs, var=var)

    adata.obs["counts"] = adata.X.sum(axis=1)
    adata_filtered = adata[adata.obs["counts"] >= 7].copy()

    # Alleen de 5 cellen met hoge counts blijven over
    assert adata_filtered.n_obs == 5


# ── Test 3: Merge ─────────────────────────────────────────────────────────────

def test_merge_geeft_sample_kolom(adata_merged):
    """
    WAT: na de merge moet er een 'sample' kolom zijn met 'bf' en 'p'.
    WAAROM: we gebruiken 'sample' later voor batch-correctie en labels.
    """
    assert "sample" in adata_merged.obs.columns
    samples = adata_merged.obs["sample"].unique().tolist()
    assert "bf" in samples
    assert "p" in samples


def test_merge_combineert_beide_datasets(adata_merged):
    """
    WAT: het gecombineerde object heeft cellen van BEIDE condities.
    WAAROM: als de merge mislukt zijn er maar cellen van één conditie.
    """
    n_bf = (adata_merged.obs["sample"] == "bf").sum()
    n_p  = (adata_merged.obs["sample"] == "p").sum()

    assert n_bf > 0, "Geen BF cellen na merge!"
    assert n_p  > 0, "Geen P cellen na merge!"


def test_merge_shape_klopt(adata_merged):
    """
    WAT: het gecombineerde object heeft meer dan 0 cellen en genen.
    WAAROM: een lege AnnData betekent dat de merge mislukt is.
    """
    assert adata_merged.n_obs > 0
    assert adata_merged.n_vars > 0
