#!/usr/bin/env python
"""Stage 3 - spatial organisation of the typed tissue.

Produces
--------
data/processed/{sample}_spatial.h5ad
results/tables/neighborhood_enrichment.csv
results/tables/niche_composition.csv
results/tables/niche_celltype_counts.csv
results/tables/immune_distance_stats.csv
results/tables/distance_profile.csv
results/figures/fig07_neighborhood_enrichment.png
results/figures/fig08_niches.png
results/figures/fig09_tumour_immune_interface.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from xtme import spatial  # noqa: E402
from xtme.config import load_config  # noqa: E402
from xtme.plotting import (  # noqa: E402
    ALARM,
    add_scalebar,
    celltype_palette,
    panel_letter,
    set_style,
    spatial_scatter,
)

TUMOUR_LABEL = "Cancer Epithelial"
IMMUNE_LABELS = ("T-cells", "B-cells", "Myeloid", "Plasmablasts")


def main() -> None:
    cfg = load_config()
    set_style()
    seed = cfg["seed"]
    sample_id = cfg["primary_sample"]
    processed, figdir, tabdir = cfg.path("processed"), cfg.path("figures"), cfg.path("tables")
    sp_cfg = cfg["spatial"]

    adata = ad.read_h5ad(processed / f"{sample_id}_annotated.h5ad")
    # Drop categories with no cells: they survive from the palette ordering and
    # would otherwise produce all-NaN rows in every pairwise statistic.
    adata.obs["cell_type"] = (adata.obs["cell_type"].astype("category")
                              .cat.remove_unused_categories())
    palette = celltype_palette(adata.obs["cell_type"])
    xy = adata.obsm["spatial"]

    spatial.build_graph(adata, method=sp_cfg["graph"], radius=sp_cfg["max_edge_um"])

    # ---- 1. neighbourhood enrichment -------------------------------------
    z = spatial.neighborhood_enrichment(adata, "cell_type", n_perms=sp_cfg["n_perms"], seed=seed)
    z.to_csv(tabdir / "neighborhood_enrichment.csv")

    # ---- 2. niches --------------------------------------------------------
    niche, profile = spatial.compute_niches(
        adata, "cell_type", k=sp_cfg["niche_k"], n_niches=sp_cfg["n_niches"], seed=seed
    )
    adata.obs["niche"] = niche.values
    profile.to_csv(tabdir / "niche_composition.csv")
    niche_names = spatial.name_niches(profile)
    pd.crosstab(adata.obs["niche"], adata.obs["cell_type"]).to_csv(
        tabdir / "niche_celltype_counts.csv")

    # ---- 3. tumour-immune interface --------------------------------------
    tumour_mask = adata.obs["cell_type"].astype(str).eq(TUMOUR_LABEL).to_numpy()
    adata.obs["dist_to_tumour"] = spatial.distance_to_compartment(adata, tumour_mask)
    stats = spatial.interface_enrichment(
        adata, "cell_type", TUMOUR_LABEL, distance_key="dist_to_tumour",
        interface_um=sp_cfg["interface_um"], n_perms=sp_cfg["n_perms"], seed=seed,
    )
    stats.to_csv(tabdir / "immune_distance_stats.csv", index=False)
    prof = spatial.distance_profile(adata, "cell_type", "dist_to_tumour")
    prof.to_csv(tabdir / "distance_profile.csv")
    print(stats.to_string(index=False))

    adata.write_h5ad(processed / f"{sample_id}_spatial.h5ad")

    # ---- Figure 7: neighbourhood enrichment ------------------------------
    order = (
        adata.obs["cell_type"].value_counts().index.tolist()
    )
    zz = z.loc[order, order]
    vmax = float(np.nanpercentile(np.abs(zz.to_numpy()), 99))
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.8),
                             gridspec_kw={"width_ratios": [1.05, 1]})
    ax = axes[0]
    im = ax.imshow(zz.to_numpy(), cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=40, ha="right")
    ax.set_yticks(range(len(order))); ax.set_yticklabels(order)
    for i in range(len(order)):
        for j in range(len(order)):
            v = zz.iat[i, j]
            ax.text(j, i, f"{v:.0f}", ha="center", va="center",
                    fontsize=plt.rcParams["xtick.labelsize"] - 1,
                    color="white" if abs(v) > vmax * 0.6 else "black")
    ax.set_title("Tumour and immune compartments avoid each other")
    plt.colorbar(im, ax=ax, fraction=0.045, pad=0.02,
                 label="neighbourhood z-score (0 = random)")

    ax = axes[1]
    pairs = []
    for i, a in enumerate(order):
        for j, b in enumerate(order):
            if j > i:
                pairs.append({"pair": f"{a} \u2013 {b}", "z": zz.iat[i, j]})
    pd.DataFrame(pairs).to_csv(tabdir / "neighborhood_pairs.csv", index=False)
    pdf = pd.DataFrame(pairs).sort_values("z")
    show = pd.concat([pdf.head(6), pdf.tail(6)])
    colours = [ALARM if v < 0 else "#2B6CB0" for v in show["z"]]
    ax.barh(np.arange(len(show)), show["z"], color=colours, height=0.72)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_yticks(np.arange(len(show))); ax.set_yticklabels(show["pair"])
    ax.set_xlabel("Neighbourhood z-score")
    ax.set_title("Most segregated and most interleaved pairs")
    for letter, ax in zip("ab", axes):
        ax.text(-0.3 if ax is axes[1] else -0.35, 1.06, letter, transform=ax.transAxes,
                fontweight="bold", fontsize=plt.rcParams["font.size"] + 2, va="bottom")
    fig.tight_layout(w_pad=3.0)
    fig.savefig(figdir / "fig07_neighborhood_enrichment.png")
    plt.close(fig)

    # ---- Figure 8: niches -------------------------------------------------
    fig = plt.figure(figsize=(12.5, 5.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.35, 1], wspace=0.22)
    ax0 = fig.add_subplot(gs[0, 0])
    niche_cats = list(profile.index)
    niche_pal = {n: c for n, c in zip(niche_cats, plt.get_cmap("tab10").colors)}
    niche_pal = {k: plt.matplotlib.colors.to_hex(v) for k, v in niche_pal.items()}
    spatial_scatter(ax0, xy, adata.obs["niche"].astype(str), categorical=True,
                    palette=niche_pal, s=0.4)
    ax0.set_title("Tissue niches: recurrent local cell-type compositions")
    add_scalebar(ax0, xy)
    # Legend labels carry the dominant population only; the second component is
    # read off panel (b), and the full two-part name stays in niche_composition.csv.
    handles = [plt.Line2D([], [], marker="o", ls="", color=niche_pal[n],
                          label=niche_names[n].split(" + ")[0], markersize=5)
               for n in niche_cats]
    # Below panel (a), not beside it: a right-hand legend lands on panel (b)'s
    # axis label, and the niche names are long enough to need the full width.
    ax0.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, -0.02),
               ncol=2, labelspacing=0.45, handletextpad=0.4, columnspacing=1.2)

    ax1 = fig.add_subplot(gs[0, 1])
    prof_plot = profile[order] if set(order) <= set(profile.columns) else profile
    bottom = np.zeros(len(prof_plot))
    for ct in prof_plot.columns:
        ax1.bar(np.arange(len(prof_plot)), prof_plot[ct], bottom=bottom,
                color=palette.get(ct, "#CBD5E0"), width=0.78, label=ct)
        bottom += prof_plot[ct].to_numpy()
    ax1.set_xticks(np.arange(len(prof_plot))); ax1.set_xticklabels(prof_plot.index)
    ax1.set_ylabel("Mean neighbourhood composition")
    ax1.set_ylim(0, 1)
    ax1.set_xlabel(f"Niche (k = {sp_cfg['niche_k']} nearest neighbours)")
    ax1.set_title("What each niche is made of")
    ax1.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=4,
               handlelength=1.0, columnspacing=1.0)
    for letter, ax in zip("ab", (ax0, ax1)):
        panel_letter(ax, letter)
    fig.savefig(figdir / "fig08_niches.png", bbox_inches="tight")
    plt.close(fig)

    # ---- Figure 9: tumour-immune interface -------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.5),
                             gridspec_kw={"width_ratios": [1.1, 1, 1]})
    ax = axes[0]
    d = np.clip(adata.obs["dist_to_tumour"].to_numpy(), 0, 200)
    h = spatial_scatter(ax, xy, d, categorical=False, cmap="cividis", s=0.35, vmin=0, vmax=200)
    ax.set_title("Distance from the tumour compartment")
    plt.colorbar(h, ax=ax, fraction=0.04, pad=0.02, label="\u00b5m (capped at 200)")
    add_scalebar(ax, xy)

    ax = axes[1]
    non_t = adata.obs[~adata.obs["cell_type"].astype(str).eq(TUMOUR_LABEL)]
    shown = [c for c in order if c != TUMOUR_LABEL][:6]
    for i, ct in enumerate(shown):
        vals = non_t.loc[non_t["cell_type"].astype(str) == ct, "dist_to_tumour"]
        parts = ax.violinplot([np.log10(vals + 1)], positions=[i], widths=0.8,
                              showextrema=False, showmedians=True)
        for body in parts["bodies"]:
            body.set_facecolor(palette.get(ct, "#CBD5E0")); body.set_alpha(0.85)
            body.set_edgecolor("none")
        parts["cmedians"].set_color("black"); parts["cmedians"].set_linewidth(1.2)
    ax.axhline(np.log10(sp_cfg["interface_um"] + 1), color=ALARM, lw=0.9, ls="--")
    ax.text(len(shown) - 0.45, np.log10(sp_cfg["interface_um"] + 1) + 0.04,
            f"{sp_cfg['interface_um']:g} \u00b5m interface", color=ALARM, ha="right",
            fontsize=plt.rcParams["legend.fontsize"])
    ax.set_xticks(range(len(shown))); ax.set_xticklabels(shown, rotation=35, ha="right")
    ax.set_ylabel("Distance to tumour (log\u2081\u2080 \u00b5m + 1)")
    ax.set_title("How close each population gets")

    ax = axes[2]
    s = stats[stats["cell_type"] != TUMOUR_LABEL].sort_values("log2_enrichment")
    colours = [palette.get(c, "#CBD5E0") for c in s["cell_type"]]
    ax.barh(np.arange(len(s)), s["log2_enrichment"], color=colours, height=0.7)
    ax.axvline(0, color="black", lw=0.8)
    for i, (_, row) in enumerate(s.iterrows()):
        star = "***" if row["q_bh"] < 0.001 else "**" if row["q_bh"] < 0.01 else \
               "*" if row["q_bh"] < 0.05 else "n.s."
        off = 0.02 if row["log2_enrichment"] >= 0 else -0.02
        ax.text(row["log2_enrichment"] + off, i, star, va="center",
                ha="left" if row["log2_enrichment"] >= 0 else "right",
                fontsize=plt.rcParams["legend.fontsize"])
    ax.set_yticks(np.arange(len(s))); ax.set_yticklabels(s["cell_type"])
    ax.set_xlabel(f"log\u2082 enrichment within {sp_cfg['interface_um']:g} \u00b5m of tumour")
    ax.set_title("Interface enrichment vs a permuted null")
    ax.margins(x=0.18)
    for letter, ax in zip("abc", axes):
        ax.text(-0.1, 1.06, letter, transform=ax.transAxes, fontweight="bold",
                fontsize=plt.rcParams["font.size"] + 2, va="bottom")
    fig.tight_layout(w_pad=2.4)
    fig.savefig(figdir / "fig09_tumour_immune_interface.png")
    plt.close(fig)

    print(f"[03] niches: {dict(adata.obs['niche'].value_counts())}")


if __name__ == "__main__":
    main()
