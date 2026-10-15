#!/usr/bin/env python
"""Stage 1 - quality control across both Xenium replicates.

Produces
--------
data/interim/{Rep1,Rep2}_qc.h5ad
results/tables/qc_scorecard.csv          per-sample cell + run metrics, verdict
results/tables/qc_thresholds.csv         the criteria that produced the verdict
results/tables/spatial_autocorrelation.csv
results/tables/replicate_gene_concordance.csv
results/figures/fig01_qc_panel.png
results/figures/fig02_spatial_qc.png
results/figures/fig03_replicate_concordance.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from xtme import io, qc  # noqa: E402
from xtme.config import load_config  # noqa: E402
from xtme.plotting import (  # noqa: E402
    ALARM,
    SAMPLE_COLORS,
    add_scalebar,
    panel_letter,
    set_style,
    spatial_scatter,
)


def main() -> None:
    cfg = load_config()
    set_style()
    rng_seed = cfg["seed"]
    interim = cfg.path("interim")
    figdir = cfg.path("figures")
    tabdir = cfg.path("tables")
    thresholds = qc.QCThresholds(**{
        k: v for k, v in cfg["qc"].items() if k in qc.QCThresholds.__dataclass_fields__
    })

    adatas, cell_rows, run_rows, moran_rows = {}, [], [], []
    for sample_id, sample_dir in cfg["samples"].items():
        print(f"[01_qc] reading {sample_id}")
        adata = io.read_xenium(sample_dir, sample_id=sample_id)
        qc.compute_cell_qc(adata, thresholds)
        cell_rows.append(qc.cell_qc_summary(adata))
        run_rows.append(qc.run_qc(adata, thresholds))

        grid = qc.grid_bin(adata, bin_um=cfg["qc"]["grid_bin_um"])
        grid.to_csv(tabdir / f"grid_qc_{sample_id}.csv", index=False)
        for metric in ("median_counts", "median_control_rate", "cell_density_per_mm2",
                       "frac_qc_pass"):
            row = qc.morans_i_grid(grid, metric, n_perm=199, seed=rng_seed)
            row["sample"] = sample_id
            moran_rows.append(row)

        adata.write_h5ad(interim / f"{sample_id}_qc.h5ad")
        adatas[sample_id] = adata

    # ---- scorecard --------------------------------------------------------
    cells = pd.DataFrame(cell_rows).set_index("sample")
    runs = pd.DataFrame(run_rows).set_index("sample")
    scorecard = runs.join(cells.drop(columns=[c for c in cells.columns if c in runs.columns]))
    scorecard.to_csv(tabdir / "qc_scorecard.csv")
    # grid_bin_um is not a pass/fail threshold but it is the spatial resolution
    # every artefact check is computed at, so it belongs in the same record.
    pd.Series({**thresholds.as_dict(), "grid_bin_um": float(cfg["qc"]["grid_bin_um"])}
              ).rename("value").to_csv(tabdir / "qc_thresholds.csv")
    pd.DataFrame(moran_rows).to_csv(tabdir / "spatial_autocorrelation.csv", index=False)
    print(scorecard[["n_cells", "median_transcripts_per_cell", "neg_control_probe_rate",
                     "frac_cells_passing", "run_verdict"]].to_string())

    # ---- Figure 1: per-cell QC -------------------------------------------
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 6.2))
    ax = axes[0, 0]
    for sid, a in adatas.items():
        counts = np.sort(a.obs["n_counts"].to_numpy())[::-1]
        ax.plot(np.arange(1, len(counts) + 1), counts, color=SAMPLE_COLORS[sid], lw=1.2, label=sid)
    ax.axhline(thresholds.min_transcripts, color=ALARM, lw=0.9, ls="--")
    ax.text(1.2, thresholds.min_transcripts * 1.15, f"filter: {thresholds.min_transcripts} transcripts",
            color=ALARM, fontsize=plt.rcParams["legend.fontsize"], va="bottom")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Cell rank"); ax.set_ylabel("Transcripts per cell")
    ax.set_title("Both replicates share one depth distribution")
    ax.legend(loc="lower left")

    ax = axes[0, 1]
    bins = np.logspace(0, np.log10(max(a.obs["n_counts"].max() for a in adatas.values())), 60)
    for sid, a in adatas.items():
        ax.hist(a.obs["n_counts"], bins=bins, histtype="step", lw=1.3,
                color=SAMPLE_COLORS[sid], density=True, label=sid)
    ax.axvline(thresholds.min_transcripts, color=ALARM, lw=0.9, ls="--")
    ax.set_xscale("log"); ax.set_xlabel("Transcripts per cell"); ax.set_ylabel("Density")
    ax.set_title("Transcript depth")

    ax = axes[0, 2]
    for sid, a in adatas.items():
        ax.hist(a.obs["n_genes"], bins=np.arange(0, a.obs["n_genes"].max() + 4, 4),
                histtype="step", lw=1.3, color=SAMPLE_COLORS[sid], density=True, label=sid)
    ax.axvline(thresholds.min_genes, color=ALARM, lw=0.9, ls="--")
    ax.set_xlabel("Genes detected per cell"); ax.set_ylabel("Density")
    ax.set_title("Panel breadth per cell")

    ax = axes[1, 0]
    primary = adatas[cfg["primary_sample"]]
    # log axes: drop the empty cells (they are counted in panel f instead)
    nz = primary.obs["n_counts"].to_numpy() > 0
    ax.hexbin(primary.obs["n_counts"][nz], primary.obs["control_rate"][nz] + 1e-4,
              xscale="log", yscale="log", gridsize=55, mincnt=1, cmap="Greys", linewidths=0)
    ax.axhline(thresholds.max_control_rate, color=ALARM, lw=0.9, ls="--")
    ax.set_xlabel("Transcripts per cell"); ax.set_ylabel("Control-feature fraction")
    ax.set_title(f"Background is concentrated in shallow cells ({cfg['primary_sample']})")

    ax = axes[1, 1]
    for sid, a in adatas.items():
        if "nucleus_ratio" in a.obs:
            ax.hist(a.obs["nucleus_ratio"], bins=np.linspace(0, 1.2, 61), histtype="step",
                    lw=1.3, color=SAMPLE_COLORS[sid], density=True, label=sid)
    ax.axvline(1.0, color=ALARM, lw=0.9, ls="--")
    ax.set_xlabel("Nucleus area / cell area"); ax.set_ylabel("Density")
    ax.set_title("Segmentation expansion")

    ax = axes[1, 2]
    flag_cols = primary.uns["qc_flag_columns"]
    frac = pd.DataFrame({sid: a.obs[flag_cols].mean() for sid, a in adatas.items()}) * 100
    frac.index = [c.replace("qc_fail_", "").replace("_", " ") for c in frac.index]
    frac = frac.sort_values(cfg["primary_sample"])
    y = np.arange(len(frac))
    for i, sid in enumerate(frac.columns):
        ax.barh(y + (i - 0.5) * 0.38, frac[sid], height=0.36,
                color=SAMPLE_COLORS[sid], label=sid)
    ax.set_yticks(y); ax.set_yticklabels(frac.index)
    ax.set_xlabel("Cells failing (%)")
    ax.set_title("Which criterion removes cells")
    ax.legend(loc="lower right")

    for letter, ax in zip("abcdef", axes.ravel()):
        ax.text(-0.14, 1.08, letter, transform=ax.transAxes, fontweight="bold",
                fontsize=plt.rcParams["font.size"] + 2, va="bottom")
    fig.tight_layout(w_pad=1.6, h_pad=2.0)
    fig.savefig(figdir / "fig01_qc_panel.png")
    plt.close(fig)

    # ---- Figure 2: spatial QC --------------------------------------------
    fig, axes = plt.subplots(2, 3, figsize=(11.5, 6.6))
    ctrl_vmax = float(np.percentile(
        np.concatenate([a.obs["control_rate"].to_numpy() for a in adatas.values()]), 99.5))
    for r, sid in enumerate(adatas):
        a = adatas[sid]
        xy = a.obsm["spatial"]
        h = spatial_scatter(axes[r, 0], xy, np.log10(a.obs["n_counts"] + 1),
                            categorical=False, cmap="viridis", s=0.35)
        axes[r, 0].set_title(f"{sid}: transcripts per cell (log\u2081\u2080)")
        plt.colorbar(h, ax=axes[r, 0], fraction=0.035, pad=0.02)
        add_scalebar(axes[r, 0], xy)

        h = spatial_scatter(axes[r, 1], xy, np.clip(a.obs["control_rate"], 0, ctrl_vmax),
                            categorical=False, cmap="magma", s=0.35, vmin=0, vmax=ctrl_vmax)
        axes[r, 1].set_title(f"{sid}: control-feature fraction")
        plt.colorbar(h, ax=axes[r, 1], fraction=0.035, pad=0.02,
                     label=f"capped at {ctrl_vmax:.3f}")

        grid = pd.read_csv(tabdir / f"grid_qc_{sid}.csv")
        sc = axes[r, 2].scatter(grid["x_um"], grid["y_um"], c=grid["cell_density_per_mm2"],
                                s=6, marker="s", cmap="cividis", lw=0)
        axes[r, 2].set_aspect("equal"); axes[r, 2].set_xticks([]); axes[r, 2].set_yticks([])
        axes[r, 2].invert_yaxis()
        for sp in axes[r, 2].spines.values():
            sp.set_visible(False)
        axes[r, 2].set_title(f"{sid}: cell density (100 \u00b5m bins)")
        plt.colorbar(sc, ax=axes[r, 2], fraction=0.035, pad=0.02, label="cells/mm\u00b2")

    for letter, ax in zip("abcdef", axes.ravel()):
        panel_letter(ax, letter)
    fig.suptitle("Depth and background vary with tissue structure, not with field-of-view position",
                 x=0.01, ha="left", fontsize=plt.rcParams["font.size"])
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(figdir / "fig02_spatial_qc.png")
    plt.close(fig)

    # ---- Figure 3: replicate concordance ---------------------------------
    pb = pd.concat([qc.pseudobulk(a) for a in adatas.values()], axis=1)
    cpm = pb / pb.sum() * 1e6
    log = np.log10(cpm + 1)
    sids = list(pb.columns)
    pearson = float(log[sids[0]].corr(log[sids[1]]))
    spearman = float(cpm[sids[0]].corr(cpm[sids[1]], method="spearman"))
    conc = pd.DataFrame({
        f"cpm_{sids[0]}": cpm[sids[0]], f"cpm_{sids[1]}": cpm[sids[1]],
        "log2_fold_change": np.log2((cpm[sids[1]] + 1) / (cpm[sids[0]] + 1)),
    })
    conc["abs_log2fc"] = conc["log2_fold_change"].abs()
    conc.sort_values("abs_log2fc", ascending=False).to_csv(
        tabdir / "replicate_gene_concordance.csv")

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.1))
    ax = axes[0]
    ax.scatter(log[sids[0]], log[sids[1]], s=7, color="#4A5568", lw=0, alpha=0.75)
    lims = [0, float(log.to_numpy().max()) * 1.05]
    ax.plot(lims, lims, color="#A0AEC0", lw=0.8, ls="--", zorder=0)
    top = conc.sort_values("abs_log2fc", ascending=False).head(6)
    for gene, row in top.iterrows():
        ax.annotate(gene, (np.log10(row[f"cpm_{sids[0]}"] + 1), np.log10(row[f"cpm_{sids[1]}"] + 1)),
                    textcoords="offset points", xytext=(4, 3),
                    fontsize=plt.rcParams["legend.fontsize"], style="italic", color=ALARM)
    ax.set_xlim(lims); ax.set_ylim(lims)
    ax.set_xlabel(f"{sids[0]} pseudobulk (log\u2081\u2080 CPM)")
    ax.set_ylabel(f"{sids[1]} pseudobulk (log\u2081\u2080 CPM)")
    ax.set_title(f"Panel-level expression replicates (r = {pearson:.3f}, \u03c1 = {spearman:.3f})")

    ax = axes[1]
    mean_expr = np.log10((cpm[sids[0]] + cpm[sids[1]]) / 2 + 1)
    ax.scatter(mean_expr, conc["log2_fold_change"], s=7, color="#4A5568", lw=0, alpha=0.75)
    ax.axhline(0, color="#A0AEC0", lw=0.8, ls="--")
    for thr in (-1, 1):
        ax.axhline(thr, color=ALARM, lw=0.7, ls=":")
    for gene, row in top.iterrows():
        ax.annotate(gene, (np.log10((row[f"cpm_{sids[0]}"] + row[f"cpm_{sids[1]}"]) / 2 + 1),
                           row["log2_fold_change"]),
                    textcoords="offset points", xytext=(4, 3),
                    fontsize=plt.rcParams["legend.fontsize"], style="italic", color=ALARM)
    ax.set_xlabel("Mean abundance (log\u2081\u2080 CPM)")
    ax.set_ylabel(f"log\u2082 {sids[1]} / {sids[0]}")
    ax.set_title("Discordant genes are rare and abundance-independent")
    for letter, ax in zip("ab", axes):
        ax.text(-0.12, 1.06, letter, transform=ax.transAxes, fontweight="bold",
                fontsize=plt.rcParams["font.size"] + 2, va="bottom")
    fig.tight_layout(w_pad=2.0)
    fig.savefig(figdir / "fig03_replicate_concordance.png")
    plt.close(fig)

    pd.Series({"pearson_log10cpm": pearson, "spearman_cpm": spearman,
               "n_genes": int(len(cpm)),
               "n_genes_abs_log2fc_gt1": int((conc["abs_log2fc"] > 1).sum())}
              ).rename("value").to_csv(tabdir / "replicate_summary.csv")
    print(f"[01_qc] replicate r={pearson:.4f} rho={spearman:.4f}")


if __name__ == "__main__":
    main()
