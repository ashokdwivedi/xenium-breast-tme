"""Synthetic fixtures so the test suite runs in seconds with no downloaded data."""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))


@pytest.fixture
def toy_xenium() -> ad.AnnData:
    """A 600-cell, 12-gene toy section with two spatially separated compartments.

    Compartment A (cells 0-299) sits at low x and expresses genes 0-5.
    Compartment B (cells 300-599) sits at high x and expresses genes 6-11.
    Ten cells are deliberately degenerate (near-zero counts, tiny area, high
    control background) so the QC flags have something to catch.
    """
    rng = np.random.default_rng(0)
    n_a = n_b = 300
    n_cells, n_genes = n_a + n_b, 12

    X = np.zeros((n_cells, n_genes))
    X[:n_a, :6] = rng.poisson(8, (n_a, 6))
    X[:n_a, 6:] = rng.poisson(1, (n_a, 6))
    X[n_a:, 6:] = rng.poisson(8, (n_b, 6))
    X[n_a:, :6] = rng.poisson(1, (n_b, 6))

    xy = np.column_stack([
        np.concatenate([rng.uniform(0, 200, n_a), rng.uniform(400, 600, n_b)]),
        rng.uniform(0, 600, n_cells),
    ])

    adata = ad.AnnData(X=sp.csr_matrix(X.astype(np.float32)))
    adata.var_names = [f"G{i}" for i in range(n_genes)]
    adata.obs_names = [f"cell{i}" for i in range(n_cells)]
    adata.obsm["spatial"] = xy
    adata.layers["counts"] = adata.X.copy()

    adata.obs["negative_control_probe_counts"] = rng.poisson(0.2, n_cells).astype(float)
    adata.obs["negative_control_codeword_counts"] = rng.poisson(0.05, n_cells).astype(float)
    adata.obs["xa_cell_area"] = rng.uniform(40, 300, n_cells)
    adata.obs["xa_nucleus_area"] = adata.obs["xa_cell_area"] * rng.uniform(0.2, 0.8, n_cells)
    adata.obs["sample"] = pd.Categorical(["Toy"] * n_cells)

    # degenerate cells
    bad = np.arange(10)
    counts = adata.layers["counts"].toarray()
    counts[bad] = 0
    counts[bad, 0] = 1
    adata.layers["counts"] = sp.csr_matrix(counts)
    adata.X = adata.layers["counts"].copy()
    adata.obs.loc[adata.obs_names[bad], "xa_cell_area"] = 5.0
    adata.obs.loc[adata.obs_names[bad], "xa_nucleus_area"] = 4.0
    adata.obs.loc[adata.obs_names[bad], "negative_control_probe_counts"] = 5.0
    return adata


@pytest.fixture
def toy_typed(toy_xenium) -> ad.AnnData:
    """``toy_xenium`` with a two-class cell type column matching the compartments."""
    labels = np.array(["TypeA"] * 300 + ["TypeB"] * 300)
    toy_xenium.obs["cell_type"] = pd.Categorical(labels)
    return toy_xenium
