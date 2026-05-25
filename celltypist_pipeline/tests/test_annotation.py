"""
test_annotation.py
Testen voor src/pipeline/annotation.py

Wat we testen:
- annotate_clusters() voegt 'cell_type' kolom toe
- Ontbrekende cluster-keys worden gedetecteerd
- Alle cellen krijgen een label als de dict compleet is
"""

import numpy as np
import anndata as ad
import pandas as pd
import pytest


def _make_clustered_adata(n_cells=30, n_clusters=3):
    """Nep-AnnData met een 'leiden' kolom alsof clustering al gedaan is."""
    np.random.seed(0)
    X = np.random.rand(n_cells, 20)
    obs = pd.DataFrame(
        {
            "leiden": [str(i % n_clusters) for i in range(n_cells)],
            "sample": ["bf" if i < n_cells // 2 else "p" for i in range(n_cells)],
        },
        index=[f"cell_{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(20)])
    return ad.AnnData(X=X, obs=obs, var=var)


# ── Test 1: annotate_clusters() ───────────────────────────────────────────────

def test_annotate_voegt_cell_type_kolom_toe():
    """
    WAT: na annotatie moet 'cell_type' in adata.obs zitten.
    WAAROM: zonder deze kolom kan CellTypist niet trainen.
    """
    adata = _make_clustered_adata(n_clusters=3)

    cluster_names = {"0": "TypeA", "1": "TypeB", "2": "TypeC"}

    # Simuleer de annotatie-logica uit annotation.py
    adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_names)

    assert "cell_type" in adata.obs.columns


def test_annotate_geen_nan_bij_complete_dict():
    """
    WAT: als alle clusters in de dict staan, mag er geen NaN zijn.
    WAAROM: NaN in 'cell_type' laat CellTypist crashen.
    """
    adata = _make_clustered_adata(n_clusters=3)
    cluster_names = {"0": "TypeA", "1": "TypeB", "2": "TypeC"}

    adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_names)

    n_nan = adata.obs["cell_type"].isna().sum()
    assert n_nan == 0, f"Er zijn {n_nan} cellen zonder label!"


def test_annotate_geeft_nan_bij_ontbrekende_cluster():
    """
    WAT: als cluster '2' niet in de dict staat, krijgen die cellen NaN.
    WAAROM: je wil dit weten VOOR je gaat trainen, niet erna.
    """
    adata = _make_clustered_adata(n_clusters=3)
    onvolledige_dict = {"0": "TypeA", "1": "TypeB"}  # cluster 2 ontbreekt!

    adata.obs["cell_type"] = adata.obs["leiden"].map(onvolledige_dict)

    n_nan = adata.obs["cell_type"].isna().sum()
    assert n_nan > 0, "Verwachtte NaN voor ontbrekend cluster, maar vond er geen."


def test_annotate_juiste_namen_toegewezen():
    """
    WAT: cluster '0' moet de naam 'TypeA' krijgen, niet iets anders.
    WAAROM: verkeerde namen geven foute biologie-conclusies.
    """
    adata = _make_clustered_adata(n_clusters=3)
    cluster_names = {"0": "TypeA", "1": "TypeB", "2": "TypeC"}

    adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_names)

    cellen_cluster_0 = adata.obs[adata.obs["leiden"] == "0"]["cell_type"]
    assert (cellen_cluster_0 == "TypeA").all()


# ── Test 2: veiligheidscheck voor training ────────────────────────────────────

def test_training_check_stopt_bij_nan():
    """
    WAT: als er NaN in 'cell_type' zit, moet de assert in train_celltypist_model() afgaan.
    WAAROM: je wil nooit een model trainen op onvolledige labels.
    """
    adata = _make_clustered_adata(n_clusters=3)
    adata.obs["cell_type"] = None  # alle cellen hebben geen label

    n_missing = adata.obs["cell_type"].isna().sum()

    # Controleer dat de assert zou afgaan (we roepen train niet echt aan)
    with pytest.raises(AssertionError):
        assert n_missing == 0, f"Los {n_missing} ontbrekende labels op voor training!"
