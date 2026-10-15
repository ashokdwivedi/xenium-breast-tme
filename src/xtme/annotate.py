"""Cell typing for a targeted spatial panel.

Two independent routes are run and then compared, because neither is
trustworthy alone on a 313-gene panel:

1. **Reference label transfer.** The Xenium cells and a breast-cancer
   scRNA-seq atlas (Wu et al. 2021) are projected into a shared PCA space
   restricted to the panel genes, batch-corrected with Harmony, and labels are
   carried over by distance-weighted k-NN with a per-cell confidence score and
   an explicit *unassigned* class.
2. **Marker scoring.** Leiden clusters are scored against canonical lineage
   marker sets and given the best-scoring identity.

Agreement between the two is evidence; disagreement localises exactly which
populations the panel cannot separate, which is itself a reportable result.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.neighbors import NearestNeighbors

# Canonical lineage markers for the breast tumour microenvironment. Only genes
# actually present on the panel are used; the rest are reported as missing so
# a weak annotation can be traced to a thin marker set rather than to biology.
MARKER_SETS: dict[str, list[str]] = {
    "Epithelial (tumour)": [
        "EPCAM", "KRT8", "KRT7", "KRT5", "KRT14", "KRT15", "CEACAM6", "ELF3",
        "ERBB2", "ESR1", "AR", "FOXA1", "TACSTD2", "CLDN4", "SERPINA3", "MLPH",
    ],
    "T cell": ["CD3D", "CD3E", "CD2", "IL7R", "CD8A", "CD8B", "TRAC", "FOXP3", "CTLA4", "PDCD1", "LCK"],
    "B cell": ["MS4A1", "CD79A", "CD79B", "BANK1", "CD19", "TCL1A"],
    "Plasma cell": ["MZB1", "DERL3", "JCHAIN", "TENT5C", "POU2AF1"],
    "Myeloid": ["CD14", "CD68", "ITGAX", "LYZ", "AIF1", "CSF1R", "FCGR3A", "ITGAM", "CD163", "C1QA", "LILRB4"],
    "Mast cell": ["KIT", "TPSAB1", "CPA3", "MS4A2", "GATA3"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "EGFL7", "CD34", "FLT1", "SOX17"],
    "Fibroblast (CAF)": ["LUM", "DCN", "POSTN", "FAP", "PDGFRA", "COL1A1", "COL1A2", "FBLN1", "SFRP4", "CXCL12"],
    "Perivascular (PVL)": ["ACTA2", "PDGFRB", "MYH11", "RGS5", "TAGLN", "NOTCH3", "MCAM"],
}


def normalise(adata: ad.AnnData, target_sum: float | None = None) -> ad.AnnData:
    """Counts-per-cell normalisation to the median depth, then ``log1p``.

    Normalising to the *median* rather than a round constant keeps the values
    on the scale of the data, which matters when a 313-gene panel gives
    per-cell totals one to two orders of magnitude below scRNA-seq.
    """
    adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    adata.layers["lognorm"] = adata.X.copy()
    return adata


def cluster(adata: ad.AnnData, n_pcs: int = 30, n_neighbors: int = 15,
            resolutions: tuple[float, ...] = (0.3, 0.5, 0.8, 1.0),
            chosen: float = 0.8, seed: int = 0) -> ad.AnnData:
    """Scale, PCA, neighbour graph, Leiden at several resolutions, UMAP.

    Every resolution is retained in ``.obs`` so the choice of ``chosen`` can be
    audited rather than asserted.
    """
    adata.X = adata.layers["lognorm"].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack", random_state=seed)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=seed)
    for res in resolutions:
        sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res}",
                     flavor="igraph", n_iterations=2, directed=False, random_state=seed)
    adata.obs["leiden"] = adata.obs[f"leiden_{chosen}"]
    sc.tl.umap(adata, random_state=seed)
    adata.X = adata.layers["lognorm"].copy()
    adata.uns["cluster_params"] = {
        "n_pcs": n_pcs, "n_neighbors": n_neighbors,
        "resolutions": list(resolutions), "chosen_resolution": chosen, "seed": seed,
    }
    return adata


def score_markers(adata: ad.AnnData, marker_sets: dict[str, list[str]] | None = None,
                  cluster_key: str = "leiden") -> pd.DataFrame:
    """Mean z-scored marker-set expression per cluster.

    Returns a cluster x lineage score matrix; ``adata.uns['marker_coverage']``
    records which markers were present on the panel.
    """
    marker_sets = marker_sets or MARKER_SETS
    present, coverage = {}, {}
    for lineage, genes in marker_sets.items():
        found = [g for g in genes if g in adata.var_names]
        coverage[lineage] = {"n_requested": len(genes), "n_found": len(found),
                             "found": found, "missing": sorted(set(genes) - set(found))}
        if found:
            present[lineage] = found
    adata.uns["marker_coverage"] = coverage

    expr = pd.DataFrame(
        np.asarray(adata[:, sorted({g for gs in present.values() for g in gs})].layers["lognorm"].todense()
                   if hasattr(adata.layers["lognorm"], "todense")
                   else adata[:, sorted({g for gs in present.values() for g in gs})].layers["lognorm"]),
        index=adata.obs_names,
        columns=sorted({g for gs in present.values() for g in gs}),
    )
    z = (expr - expr.mean()) / expr.std(ddof=0).replace(0, np.nan)
    groups = adata.obs[cluster_key].astype(str)
    scores = pd.DataFrame(
        {lineage: z[genes].mean(axis=1) for lineage, genes in present.items()}
    ).groupby(groups).mean()
    scores.index.name = cluster_key
    return scores


def label_clusters_by_marker(scores: pd.DataFrame, margin: float = 0.10) -> pd.DataFrame:
    """Assign each cluster its top-scoring lineage, with an ambiguity flag.

    ``margin`` is the minimum gap between the best and second-best lineage
    score; clusters below it are kept but marked ambiguous rather than being
    silently given a confident label.
    """
    ranked = scores.apply(lambda r: r.sort_values(ascending=False), axis=1)
    best = scores.idxmax(axis=1)
    top2 = np.sort(scores.to_numpy(), axis=1)[:, -2:]
    gap = top2[:, 1] - top2[:, 0]
    return pd.DataFrame({
        "marker_label": best,
        "top_score": scores.max(axis=1).to_numpy(),
        "score_gap": gap,
        "ambiguous": gap < margin,
        "runner_up": [ranked.iloc[i].index[1] for i in range(len(ranked))],
    }, index=scores.index)


def transfer_labels(
    query: ad.AnnData,
    reference: ad.AnnData,
    label_keys: tuple[str, ...] = ("celltype_major", "celltype_minor"),
    n_pcs: int = 30,
    k: int = 25,
    min_confidence: float = 0.5,
    seed: int = 0,
) -> pd.DataFrame:
    """Harmony-integrated distance-weighted k-NN label transfer.

    Parameters
    ----------
    query, reference
        Must already share an identical, identically ordered ``var_names``
        (the panel genes) and carry a ``counts`` layer.
    k
        Reference neighbours per query cell.
    min_confidence
        Minimum summed neighbour weight for the winning label. Cells below it
        are labelled ``Unassigned`` instead of being forced into the nearest
        class, which is the honest outcome when a panel cannot resolve a cell.

    Returns
    -------
    DataFrame indexed like ``query`` with ``<key>``, ``<key>_confidence`` and
    ``<key>_entropy`` columns.
    """
    import harmonypy

    if list(query.var_names) != list(reference.var_names):
        raise ValueError("query and reference must share identical var_names")

    joint = ad.concat(
        {"reference": reference, "query": query}, label="batch",
        join="outer", index_unique="-",
    )
    joint.layers["counts"] = joint.X.copy()
    normalise(joint)
    sc.pp.scale(joint, max_value=10)
    sc.tl.pca(joint, n_comps=n_pcs, svd_solver="arpack", random_state=seed)

    harmony = harmonypy.run_harmony(
        joint.obsm["X_pca"], joint.obs, ["batch"], max_iter_harmony=20, random_state=seed
    )
    # harmonypy has returned Z_corr as (n_pcs, n_cells) in some releases and
    # (n_cells, n_pcs) in others; orient by the dimension we know.
    Z = np.asarray(harmony.Z_corr)
    emb = Z if Z.shape[0] == joint.n_obs else Z.T
    is_ref = (joint.obs["batch"] == "reference").to_numpy()
    ref_emb, qry_emb = emb[is_ref], emb[~is_ref]

    nn = NearestNeighbors(n_neighbors=k, algorithm="auto").fit(ref_emb)
    dist, idx = nn.kneighbors(qry_emb)
    # Gaussian kernel with a per-cell bandwidth = distance to the k-th neighbour
    sigma = np.maximum(dist[:, -1:], 1e-12)
    weights = np.exp(-(dist / sigma) ** 2)
    weights /= weights.sum(axis=1, keepdims=True)

    out = pd.DataFrame(index=query.obs_names)
    for key in label_keys:
        if key not in reference.obs:
            continue
        ref_labels = reference.obs[key].astype(str).to_numpy()
        classes, codes = np.unique(ref_labels, return_inverse=True)
        votes = np.zeros((len(qry_emb), len(classes)), dtype=np.float32)
        neighbour_codes = codes[idx]
        for col in range(k):
            np.add.at(votes, (np.arange(len(qry_emb)), neighbour_codes[:, col]), weights[:, col])
        winner = votes.argmax(axis=1)
        confidence = votes.max(axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            p = np.where(votes > 0, votes, 1.0)
            entropy = -(votes * np.log(p)).sum(axis=1)
        labels = classes[winner].astype(object)
        labels[confidence < min_confidence] = "Unassigned"
        out[key] = pd.Categorical(labels)
        out[f"{key}_confidence"] = confidence
        out[f"{key}_entropy"] = entropy

    out.attrs["harmony_k"] = k
    out.attrs["min_confidence"] = min_confidence
    return out


def concordance(a: pd.Series, b: pd.Series, normalize: str = "index") -> pd.DataFrame:
    """Row-normalised contingency table between two labellings of the same cells."""
    return pd.crosstab(a.astype(str), b.astype(str), normalize=normalize)
