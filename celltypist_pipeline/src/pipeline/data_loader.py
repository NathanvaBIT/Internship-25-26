"""
data_loader.py
Laadt de BF en P h5ad bestanden en geeft een gecombineerd AnnData object terug.
"""

import scanpy as sc
import anndata as ad


def load_and_merge(bf_path: str, p_path: str) -> ad.AnnData:
    """
    Laad BF en P datasets, filter ribosomale genen, pas count-filter toe,
    normaliseer en merge tot één AnnData object.

    Parameters
    ----------
    bf_path : pad naar BF_count_mat.h5ad
    p_path  : pad naar P_count_mat.h5ad

    Returns
    -------
    adata_bacteria : gecombineerd AnnData object met 'sample' kolom (bf / p)
    """
    # --- Laden ---
    adata_bf = sc.read_h5ad(bf_path)
    adata_p  = sc.read_h5ad(p_path)

    # --- Ribo-filter ---
    ribo_prefixes = ("rps", "rpl", "rpm", "rrs", "rrl", "rrf")

    adata_bf.var_names_make_unique()
    adata_bf.var["ribo"] = adata_bf.var_names.str.lower().str.startswith(ribo_prefixes)
    adata_bf = adata_bf[:, ~adata_bf.var["ribo"]].copy()

    adata_p.var_names_make_unique()
    adata_p.var["ribo"] = adata_p.var_names.str.lower().str.startswith(ribo_prefixes)
    adata_p = adata_p[:, ~adata_p.var["ribo"]].copy()

    # --- Count-filter (cutoffs uit het paper) ---
    adata_bf.obs["counts"] = adata_bf.X.sum(axis=1)
    adata_bf = adata_bf[adata_bf.obs["counts"] >= 7].copy()

    adata_p.obs["counts"] = adata_p.X.sum(axis=1)
    adata_p = adata_p[adata_p.obs["counts"] >= 28].copy()

    # --- Normalisatie ---
    sc.pp.normalize_total(adata_bf, target_sum=1e4)
    sc.pp.log1p(adata_bf)

    sc.pp.normalize_total(adata_p, target_sum=1e4)
    sc.pp.log1p(adata_p)

    # --- Samenvoegen ---
    adata_bacteria = ad.concat(
        {"bf": adata_bf, "p": adata_p},
        label="sample",
        join="inner"
    )
    adata_bacteria.obs_names_make_unique()

    print(f"Dataset geladen — shape: {adata_bacteria.shape}")
    print(f"Labels aanwezig: {adata_bacteria.obs['sample'].unique().tolist()}")

    return adata_bacteria
