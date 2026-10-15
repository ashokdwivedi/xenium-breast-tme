#!/usr/bin/env python
"""Stage 2 - cluster the primary replicate and type its cells two ways.

Route A: Harmony-integrated k-NN label transfer from the Wu et al. 2021
         breast-cancer scRNA-seq atlas (GSE176078), restricted to panel genes.
Route B: canonical marker scoring of Leiden clusters.

The two are compared rather than merged, because where they disagree is where
a 313-gene panel stops resolving cell identity.

Produces
--------
data/processed/{sample}_annotated.h5ad
results/tables/celltype_counts.csv
results/tables/marker_coverage.csv
results/tables/cluster_marker_labels.csv
results/tables/celltype_concordance.csv
results/figures/fig04_umap_clusters.png
results/figures/fig05_celltypes.png
results/figures/fig06_annotation_quality.png
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from xtme import annotate, io  # noqa: E402
from xtme.config import load_config  # noqa: E402
from xtme.plotting import (  # noqa: E402
    add_scalebar,
    celltype_palette,
    panel_letter,
    set_style,
    spatial_scatter,
)


def main() -> None:
    cfg = load_config()
    set_style()
    seed = cfg["seed"]
    sample_id = cfg["primary_sample"]
    interim, processed = cfg.path("interim"), cfg.path("processed")
    figdir, tabdir = cfg.path("figures"), cfg.path("tables")

    # ---- load QC'd data and apply the cell filter ------------------------
    import anndata as ad
    adata = ad.read_h5ad(interim / f"{sample_id}_qc.h5ad")
    n_before = adata.n_obs
    adata = adata[adata.obs["qc_pass"].to_numpy()].copy()
    print(f"[02] {sample_id}: {adata.n_obs}/{n_before} cells pass QC")

    annotate.normalise(adata)
    annotate.cluster(
        adata, n_pcs=cfg["cluster"]["n_pcs"], n_neighbors=cfg["cluster"]["n_neighbors"],
        resolutions=tuple(cfg["cluster"]["resolutions"]),
        chosen=cfg["cluster"]["chosen_resolution"], seed=seed,
    )
    print(f"[02] leiden clusters: {adata.obs['leiden'].nunique()}")

    # ---- Route B: marker scoring -----------------------------------------
    scores = annotate.score_markers(adata)
    labels = annotate.label_clusters_by_marker(scores)
    scores.to_csv(tabdir / "cluster_marker_scores.csv")
    labels.to_csv(tabdir / "cluster_marker_labels.csv")
    pd.DataFrame([
        {"lineage": k, "n_requested": v["n_requested"], "n_found": v["n_found"],
         "found": ";".join(v["found"]), "missing": ";".join(v["missing"])}
        for k, v in adata.uns["marker_coverage"].items()
    ]).to_csv(tabdir / "marker_coverage.csv", index=False)
    adata.obs["marker_label"] = adata.obs["leiden"].astype(str).map(labels["marker_label"]).astype("category")
    adata.obs["marker_ambiguous"] = adata.obs["leiden"].astype(str).map(labels["ambiguous"]).astype(bool)

    # ---- Route A: reference label transfer -------------------------------
    panel = list(adata.var_names)
    ref = io.read_reference_subset(cfg["reference"]["dir"], keep_genes=panel)
    print(f"[02] reference: {ref.n_obs} cells x {ref.n_vars} panel genes")

    # Down-sample abundant reference labels so the k-NN vote is not decided by
    # class prior alone, and so Harmony stays within memory.
    cap = cfg["reference"]["max_cells_per_label"]
    key0 = cfg["reference"]["label_keys"][0]
    ref = ref[ref.obs[key0].notna()].copy()
    rng = np.random.default_rng(seed)
    keep = np.concatenate([
        rng.choice(np.flatnonzero((ref.obs[key0] == lab).to_numpy()),
                   size=min(cap, int((ref.obs[key0] == lab).sum())), replace=False)
        for lab in ref.obs[key0].astype(str).unique()
    ])
    ref = ref[np.sort(keep)].copy()
    print(f"[02] reference after capping at {cap}/label: {ref.n_obs} cells")

    transferred = annotate.transfer_labels(
        adata, ref, label_keys=tuple(cfg["reference"]["label_keys"]),
        n_pcs=cfg["cluster"]["n_pcs"], k=cfg["transfer"]["k"],
        min_confidence=cfg["transfer"]["min_confidence"], seed=seed,
    )
    for col in transferred.columns:
        adata.obs[col] = transferred[col].values

    # ---- consensus: majority transferred label per Leiden cluster --------
    major = cfg["reference"]["label_keys"][0]
    cluster_major = (
        adata.obs.groupby("leiden", observed=True)[major]
        .agg(lambda s: s.value_counts().idxmax())
    )
    adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_major).astype("category")

    counts = (
        adata.obs.groupby("cell_type", observed=True)
        .size().rename("n_cells").to_frame()
        .assign(pct=lambda d: d["n_cells"] / d["n_cells"].sum() * 100)
        .sort_values("n_cells", ascending=False)
    )
    counts.to_csv(tabdir / "celltype_counts.csv")
    print(counts.to_string())

    conc = annotate.concordance(adata.obs["cell_type"], adata.obs["marker_label"])
    conc.to_csv(tabdir / "celltype_concordance.csv")

    cluster_table = pd.DataFrame({
        "n_cells": adata.obs.groupby("leiden", observed=True).size(),
        "transferred_majority": cluster_major,
        "marker_label": labels["marker_label"],
        "marker_ambiguous": labels["ambiguous"],
        "mean_transfer_confidence": adata.obs.groupby("leiden", observed=True)[f"{major}_confidence"].mean(),
        "median_n_counts": adata.obs.groupby("leiden", observed=True)["n_counts"].median(),
    })
    cluster_table["routes_agree"] = [
        _same(a, b) for a, b in zip(cluster_table["transferred_majority"], cluster_table["marker_label"])
    ]
    cluster_table.to_csv(tabdir / "cluster_annotation_table.csv")

    adata.write_h5ad(processed / f"{sample_id}_annotated.h5ad")

    # ---- figures ---------------------------------------------------------
    palette = celltype_palette(adata.obs["cell_type"])
    xy = adata.obsm["spatial"]

    # Fig 4: Leiden clusters, UMAP + tissue
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 5.0))
    sc.pl.umap(adata, color="leiden", ax=axes[0], show=False, legend_loc="on data",
               legend_fontsize=6, size=2, frameon=False, title="")
    axes[0].set_title(f"Leiden (res {cfg['cluster']['chosen_resolution']}) separates "
                      f"{adata.obs['leiden'].nunique()} transcriptional states")
    lei_pal = {c: col for c, col in zip(adata.obs["leiden"].cat.categories,
                                        adata.uns["leiden_colors"])}
    spatial_scatter(axes[1], xy, adata.obs["leiden"], categorical=True, palette=lei_pal, s=0.35)
    axes[1].set_title("The same clusters, placed back in tissue")
    add_scalebar(axes[1], xy)
    for letter, ax in zip("ab", axes):
        panel_letter(ax, letter)
    fig.tight_layout(w_pad=1.5)
    fig.savefig(figdir / "fig04_umap_clusters.png")
    plt.close(fig)

    # Fig 5: cell types
    fig = plt.figure(figsize=(12.0, 5.4))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1.35, 0.8], wspace=0.30)
    ax0 = fig.add_subplot(gs[0, 0])
    for lab, colour in palette.items():
        m = (adata.obs["cell_type"].astype(str) == lab).to_numpy()
        ax0.scatter(adata.obsm["X_umap"][m, 0], adata.obsm["X_umap"][m, 1],
                    s=1.2, c=colour, lw=0, rasterized=True)
    ax0.set_xticks([]); ax0.set_yticks([])
    for sp in ax0.spines.values():
        sp.set_visible(False)
    ax0.set_title("UMAP")
    _umap_axis_cue(ax0)

    ax1 = fig.add_subplot(gs[0, 1])
    spatial_scatter(ax1, xy, adata.obs["cell_type"], categorical=True, palette=palette, s=0.4)
    ax1.set_title("Cell types in tissue")
    add_scalebar(ax1, xy)
    # No legend: panel (c) lists every cell type against the same colours, so a
    # legend here would duplicate it and crowd the tissue map.
    ax1.get_legend().remove() if ax1.get_legend() else None

    ax2 = fig.add_subplot(gs[0, 2])
    order = counts.index.tolist()[::-1]
    ax2.barh(np.arange(len(order)), counts.loc[order, "pct"],
             color=[palette[o] for o in order], height=0.72)
    for i, o in enumerate(order):
        ax2.text(counts.loc[o, "pct"] + 0.6, i, f"{counts.loc[o, 'pct']:.1f}%",
                 va="center", fontsize=plt.rcParams["legend.fontsize"])
    ax2.set_yticks(np.arange(len(order))); ax2.set_yticklabels(order)
    ax2.set_xlabel("Cells (%)")
    ax2.set_xlim(0, counts["pct"].max() * 1.25)
    ax2.set_title("Composition")
    for letter, ax in zip("abc", (ax0, ax1, ax2)):
        panel_letter(ax, letter)
    fig.suptitle(f"{sample_id}: {adata.n_obs:,} QC-passing cells typed against a breast-cancer "
                 "scRNA-seq reference", x=0.01, ha="left", fontsize=plt.rcParams["font.size"])
    fig.savefig(figdir / "fig05_celltypes.png", bbox_inches="tight")
    plt.close(fig)

    # Fig 6: annotation quality
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.3))
    ax = axes[0]
    conf = adata.obs[f"{major}_confidence"]
    for lab in counts.index:
        m = (adata.obs["cell_type"].astype(str) == lab).to_numpy()
        ax.scatter(np.random.default_rng(0).normal(list(counts.index).index(lab), 0.09, m.sum()),
                   conf[m], s=0.6, c=palette[lab], lw=0, alpha=0.25, rasterized=True)
        ax.plot([list(counts.index).index(lab) - 0.3, list(counts.index).index(lab) + 0.3],
                [np.median(conf[m])] * 2, color="black", lw=1.4)
    ax.axhline(cfg["transfer"]["min_confidence"], color="#D62728", lw=0.9, ls="--")
    ax.set_xticks(np.arange(len(counts.index)))
    ax.set_xticklabels(counts.index, rotation=35, ha="right")
    ax.set_ylabel("k-NN label confidence")
    ax.set_title("Transfer confidence by assigned type")

    ax = axes[1]
    im = ax.imshow(conc.to_numpy(), cmap="Blues", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(conc.shape[1])); ax.set_xticklabels(conc.columns, rotation=35, ha="right")
    ax.set_yticks(range(conc.shape[0])); ax.set_yticklabels(conc.index)
    for i in range(conc.shape[0]):
        for j in range(conc.shape[1]):
            v = conc.iat[i, j]
            if v >= 0.01:
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=plt.rcParams["xtick.labelsize"] - 1,
                        color="white" if v > 0.55 else "black")
    ax.set_xlabel("Marker-based label"); ax.set_ylabel("Reference-transferred label")
    ax.set_title("Where the two routes agree")
    plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02, label="row fraction")

    ax = axes[2]
    cov = pd.read_csv(tabdir / "marker_coverage.csv")
    cov = cov.sort_values("n_found")
    ax.barh(np.arange(len(cov)), cov["n_found"], color="#4A5568", height=0.7)
    ax.barh(np.arange(len(cov)), cov["n_requested"] - cov["n_found"],
            left=cov["n_found"], color="#CBD5E0", height=0.7)
    ax.set_yticks(np.arange(len(cov))); ax.set_yticklabels(cov["lineage"])
    ax.set_xlabel("Canonical markers")
    ax.set_title("Panel coverage limits what is resolvable")
    ax.text(0.98, 0.04, "dark = on panel   light = absent", transform=ax.transAxes,
            ha="right", fontsize=plt.rcParams["legend.fontsize"], color="#4A5568")
    for letter, ax in zip("abc", axes):
        ax.text(-0.14, 1.06, letter, transform=ax.transAxes, fontweight="bold",
                fontsize=plt.rcParams["font.size"] + 2, va="bottom")
    fig.tight_layout(w_pad=2.2)
    fig.savefig(figdir / "fig06_annotation_quality.png")
    plt.close(fig)

    summary = {
        "sample": sample_id,
        "n_cells_qc_pass": int(adata.n_obs),
        "n_leiden_clusters": int(adata.obs["leiden"].nunique()),
        "n_reference_cells": int(ref.n_obs),
        "frac_unassigned": float((adata.obs[major].astype(str) == "Unassigned").mean()),
        "median_transfer_confidence": float(conf.median()),
        "frac_clusters_routes_agree": float(cluster_table["routes_agree"].mean()),
    }
    (tabdir / "annotation_summary.json").write_text(json.dumps(summary, indent=2))
    print("[02]", summary)


def _same(transferred: str, marker: str) -> bool:
    """Map the reference vocabulary onto the marker vocabulary for comparison."""
    alias = {
        "Cancer Epithelial": "Epithelial (tumour)", "Normal Epithelial": "Epithelial (tumour)",
        "T-cells": "T cell", "B-cells": "B cell", "Plasmablasts": "Plasma cell",
        "CAFs": "Fibroblast (CAF)", "PVL": "Perivascular (PVL)",
        "Myeloid": "Myeloid", "Endothelial": "Endothelial",
    }
    return alias.get(str(transferred), str(transferred)) == str(marker)


def _umap_axis_cue(ax) -> None:
    ax.annotate("", xy=(0.14, 0.02), xytext=(0.02, 0.02), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="->", lw=0.8, color="#4A5568"))
    ax.annotate("", xy=(0.02, 0.14), xytext=(0.02, 0.02), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="->", lw=0.8, color="#4A5568"))
    ax.text(0.16, 0.02, "UMAP1", transform=ax.transAxes, fontsize=plt.rcParams["xtick.labelsize"],
            color="#4A5568", va="center")
    ax.text(0.02, 0.16, "UMAP2", transform=ax.transAxes, fontsize=plt.rcParams["xtick.labelsize"],
            color="#4A5568", ha="left", va="bottom")


if __name__ == "__main__":
    main()
