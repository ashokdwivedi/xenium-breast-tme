"""Spatial organisation: neighbourhood enrichment, tissue niches, interfaces.

The unit of analysis here is not the cell but the *company a cell keeps*. Three
complementary views:

* :func:`neighborhood_enrichment` - which cell-type pairs touch more (or less)
  often than a permuted null allows.
* :func:`compute_niches` - recurrent local compositions, found by clustering
  each cell's neighbourhood composition vector. These are the tissue
  "ecosystems": tumour nest, stromal band, immune aggregate, vasculature.
* :func:`distance_to_compartment` / :func:`interface_enrichment` - how far each
  cell sits from a reference compartment (here, tumour), and whether a subset
  is enriched at the interface relative to a label-permuted null.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import squidpy as sq
from sklearn.cluster import KMeans
from sklearn.neighbors import NearestNeighbors


def build_graph(adata: ad.AnnData, method: str = "delaunay", n_neighs: int = 10,
                radius: float | None = None) -> ad.AnnData:
    """Build the spatial neighbour graph used by every downstream statistic.

    ``delaunay`` adapts to local cell density, which is the right default in
    tissue where tumour nests are far denser than stroma; a fixed-``k`` graph
    silently spans different physical distances in the two compartments.
    """
    if method == "delaunay":
        sq.gr.spatial_neighbors(adata, delaunay=True, coord_type="generic")
        # Delaunay edges can span large gaps at the tissue border; cap them at
        # a generous physical distance so holes do not create false adjacency.
        if radius is not None:
            G, D = adata.obsp["spatial_connectivities"], adata.obsp["spatial_distances"]
            mask = (D > radius)
            G = G.tolil(); D = D.tolil()
            G[mask] = 0; D[mask] = 0
            adata.obsp["spatial_connectivities"] = G.tocsr()
            adata.obsp["spatial_distances"] = D.tocsr()
            adata.obsp["spatial_connectivities"].eliminate_zeros()
    elif method == "knn":
        sq.gr.spatial_neighbors(adata, n_neighs=n_neighs, coord_type="generic")
    elif method == "radius":
        sq.gr.spatial_neighbors(adata, radius=radius, coord_type="generic")
    else:
        raise ValueError(f"unknown graph method: {method}")
    adata.uns["spatial_graph_params"] = {"method": method, "n_neighs": n_neighs, "radius": radius}
    return adata


def neighborhood_enrichment(adata: ad.AnnData, cluster_key: str,
                            n_perms: int = 1000, seed: int = 0) -> pd.DataFrame:
    """Permutation z-scores for cell-type co-occurrence on the spatial graph."""
    # n_jobs=1: squidpy's multiprocessing backend opens a unix socket in TMPDIR,
    # which fails on long temp paths, and the permutation loop is numba-compiled
    # so the single-process path is not the bottleneck.
    sq.gr.nhood_enrichment(adata, cluster_key=cluster_key, n_perms=n_perms, seed=seed,
                           n_jobs=1, show_progress_bar=False)
    cats = list(adata.obs[cluster_key].cat.categories)
    z = adata.uns[f"{cluster_key}_nhood_enrichment"]["zscore"]
    return pd.DataFrame(z, index=cats, columns=cats)


def neighbourhood_composition(adata: ad.AnnData, cluster_key: str,
                              k: int = 20) -> pd.DataFrame:
    """Fraction of each cell type among each cell's ``k`` nearest neighbours.

    Uses a fixed-``k`` spatial neighbourhood (not the Delaunay graph) so every
    composition vector is built from the same number of cells and the
    downstream clustering is not confounded by local density.
    """
    xy = adata.obsm["spatial"]
    nn = NearestNeighbors(n_neighbors=k + 1).fit(xy)
    _, idx = nn.kneighbors(xy)
    idx = idx[:, 1:]  # drop self

    labels = adata.obs[cluster_key].astype("category")
    cats = list(labels.cat.categories)
    codes = labels.cat.codes.to_numpy()

    comp = np.zeros((adata.n_obs, len(cats)), dtype=np.float32)
    neighbour_codes = codes[idx]
    for col in range(k):
        np.add.at(comp, (np.arange(adata.n_obs), neighbour_codes[:, col]), 1.0)
    comp /= k
    return pd.DataFrame(comp, index=adata.obs_names, columns=cats)


def compute_niches(adata: ad.AnnData, cluster_key: str, k: int = 20,
                   n_niches: int = 8, seed: int = 0) -> tuple[pd.Series, pd.DataFrame]:
    """Cluster neighbourhood composition vectors into recurrent tissue niches.

    Returns the per-cell niche assignment and the niche x cell-type mean
    composition matrix that defines what each niche *is*.
    """
    comp = neighbourhood_composition(adata, cluster_key=cluster_key, k=k)
    km = KMeans(n_clusters=n_niches, random_state=seed, n_init=10).fit(comp.to_numpy())
    niche = pd.Series(
        pd.Categorical([f"N{i}" for i in km.labels_],
                       categories=[f"N{i}" for i in range(n_niches)]),
        index=adata.obs_names, name="niche",
    )
    profile = comp.groupby(niche.values).mean()
    profile.index.name = "niche"
    adata.obsm["nhood_composition"] = comp.to_numpy()
    adata.uns["niche_params"] = {"k": k, "n_niches": n_niches, "seed": seed,
                                 "cluster_key": cluster_key}
    return niche, profile


def name_niches(profile: pd.DataFrame, top_n: int = 2) -> dict[str, str]:
    """Human-readable niche names from the top contributing cell types."""
    names = {}
    for niche, row in profile.iterrows():
        top = row.sort_values(ascending=False).head(top_n)
        names[str(niche)] = f"{niche}: " + " + ".join(
            f"{ct} ({frac:.0%})" for ct, frac in top.items()
        )
    return names


# ---------------------------------------------------------------------------
# tumour-immune interface
# ---------------------------------------------------------------------------

def distance_to_compartment(adata: ad.AnnData, mask: np.ndarray) -> np.ndarray:
    """Euclidean distance (um) from every cell to the nearest cell in ``mask``.

    Cells inside the compartment get distance 0 by construction (they are their
    own nearest member), so downstream tests should exclude them explicitly.
    """
    if mask.sum() == 0:
        raise ValueError("compartment mask is empty")
    xy = adata.obsm["spatial"]
    nn = NearestNeighbors(n_neighbors=1).fit(xy[mask])
    dist, _ = nn.kneighbors(xy)
    return dist.ravel()


def interface_enrichment(adata: ad.AnnData, celltype_key: str, target: str,
                         distance_key: str = "dist_to_tumour",
                         interface_um: float = 30.0, n_perms: int = 1000,
                         seed: int = 0) -> pd.DataFrame:
    """Test whether each non-target cell type is enriched at the interface.

    For each cell type, the statistic is the fraction of its cells lying within
    ``interface_um`` of the target compartment. The null permutes cell-type
    labels among all non-target cells, holding positions fixed, so it controls
    for the overall geometry of the tissue and for how much of the section is
    near tumour at all. Reports the observed fraction, the null mean, a
    log2 enrichment and a two-sided permutation p-value.
    """
    rng = np.random.default_rng(seed)
    obs = adata.obs
    is_target = obs[celltype_key].astype(str).eq(target).to_numpy()
    pool = ~is_target
    labels = obs.loc[pool, celltype_key].astype(str).to_numpy()
    near = (obs.loc[pool, distance_key].to_numpy() <= interface_um)

    rows = []
    cats = [c for c in pd.unique(labels)]
    observed = {c: near[labels == c].mean() for c in cats}

    null = {c: np.empty(n_perms) for c in cats}
    for p in range(n_perms):
        shuffled = rng.permutation(labels)
        for c in cats:
            null[c][p] = near[shuffled == c].mean()

    for c in cats:
        o, nl = observed[c], null[c]
        mu = float(nl.mean())
        pval = (np.sum(np.abs(nl - mu) >= abs(o - mu)) + 1) / (n_perms + 1)
        rows.append({
            "cell_type": c,
            "n_cells": int((labels == c).sum()),
            "frac_within_interface": float(o),
            "null_mean": mu,
            "null_sd": float(nl.std(ddof=1)),
            "log2_enrichment": float(np.log2((o + 1e-9) / (mu + 1e-9))),
            "z_score": float((o - mu) / nl.std(ddof=1)) if nl.std(ddof=1) > 0 else np.nan,
            "p_perm": float(pval),
        })
    out = pd.DataFrame(rows).sort_values("log2_enrichment", ascending=False)
    out["interface_um"] = interface_um
    # Benjamini-Hochberg across cell types
    m = len(out)
    order = out["p_perm"].rank(method="first").to_numpy()
    out["q_bh"] = np.minimum(1.0, out["p_perm"].to_numpy() * m / order)
    return out.reset_index(drop=True)


def distance_profile(adata: ad.AnnData, celltype_key: str, distance_key: str,
                     bins: np.ndarray | None = None) -> pd.DataFrame:
    """Composition of the tissue as a function of distance from the compartment."""
    bins = bins if bins is not None else np.array([0, 10, 20, 30, 50, 75, 100, 150, 200, 300, 1e6])
    d = adata.obs[distance_key].to_numpy()
    binned = pd.cut(d, bins=bins, right=False)
    tab = pd.crosstab(binned, adata.obs[celltype_key].astype(str), normalize="index")
    tab.index.name = f"{distance_key}_bin"
    return tab
