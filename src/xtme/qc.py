"""Quality control for Xenium in situ runs.

Three levels, deliberately kept separate because they fail for different
reasons and have different remedies:

``cell_qc``
    Per-cell metrics. Drives cell filtering.
``run_qc``
    Per-sample/run summaries, including the negative-control background rate
    that distinguishes a decoding problem from a tissue problem.
``spatial_qc``
    Position-aware checks. A run can look healthy in every histogram and still
    have a folded region, a detached tissue edge or a dead field of view. These
    are only visible once metrics are binned back onto tissue coordinates.

Thresholds live in :class:`QCThresholds` so they are explicit, versioned and
overridable per study rather than scattered through a notebook.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


@dataclass(frozen=True)
class QCThresholds:
    """Cell- and run-level acceptance criteria.

    Defaults are chosen for a 313-plex FFPE panel and are justified in
    ``docs/pages/qc.md``. They are starting points for a panel of this size,
    not universal constants: a 5000-plex panel shifts every count-based
    threshold upward, and a sparse tissue shifts the area thresholds.
    """

    # --- cell level ---
    min_transcripts: int = 10        # below this a cell cannot be typed reliably
    min_genes: int = 5               # guards against one gene carrying a cell
    min_cell_area: float = 15.0      # um^2; smaller objects are segmentation debris
    max_cell_area: float = 1500.0    # um^2; larger objects are usually merged cells
    max_control_rate: float = 0.05   # per-cell fraction of counts from controls
    max_nucleus_ratio: float = 1.0   # nucleus cannot exceed its own cell

    # --- run level (pass/fail on the whole sample) ---
    max_run_neg_probe_rate: float = 0.01      # 1% of panel counts
    max_run_neg_codeword_rate: float = 0.005  # 0.5% of panel counts
    min_median_transcripts: float = 20.0
    min_cells_passing: float = 0.70           # fraction surviving cell filters

    def as_dict(self) -> dict:
        return asdict(self)


CONTROL_COUNT_COLS = (
    "negative_control_probe_counts",
    "negative_control_codeword_counts",
    "unassigned_codeword_counts",
    "genomic_control_counts",
    "deprecated_codeword_counts",
)


def compute_cell_qc(adata: ad.AnnData, thresholds: QCThresholds | None = None) -> ad.AnnData:
    """Add per-cell QC columns and boolean flags to ``adata.obs`` (in place).

    Columns added
    -------------
    ``n_counts``, ``n_genes``       panel transcripts and distinct genes per cell
    ``control_counts``              summed counts across all control classes
    ``control_rate``                control / (control + panel)
    ``nucleus_ratio``               nucleus_area / cell_area
    ``transcript_density``          panel transcripts per um^2 of segmented cell
    ``qc_fail_*``                   one boolean per criterion
    ``qc_pass``                     conjunction of all criteria
    """
    thresholds = thresholds or QCThresholds()
    X = adata.layers.get("counts", adata.X)
    X = sp.csr_matrix(X) if not sp.issparse(X) else X

    adata.obs["n_counts"] = np.asarray(X.sum(axis=1)).ravel()
    adata.obs["n_genes"] = np.asarray((X > 0).sum(axis=1)).ravel()

    present = [c for c in CONTROL_COUNT_COLS if c in adata.obs]
    adata.obs["control_counts"] = (
        adata.obs[present].sum(axis=1).to_numpy() if present else np.zeros(adata.n_obs)
    )
    denom = adata.obs["n_counts"].to_numpy() + adata.obs["control_counts"].to_numpy()
    adata.obs["control_rate"] = np.divide(
        adata.obs["control_counts"].to_numpy(), denom,
        out=np.zeros(adata.n_obs), where=denom > 0,
    )

    area = adata.obs.get("xa_cell_area")
    nuc = adata.obs.get("xa_nucleus_area")
    if area is not None:
        area = area.to_numpy(dtype=float)
        adata.obs["transcript_density"] = np.divide(
            adata.obs["n_counts"].to_numpy(), area,
            out=np.zeros(adata.n_obs), where=area > 0,
        )
        if nuc is not None:
            nuc = nuc.to_numpy(dtype=float)
            adata.obs["nucleus_ratio"] = np.divide(
                nuc, area, out=np.zeros(adata.n_obs), where=area > 0
            )
            adata.obs["no_nucleus"] = nuc <= 0

    # --- flags -----------------------------------------------------------
    t = thresholds
    flags = {
        "qc_fail_low_counts": adata.obs["n_counts"] < t.min_transcripts,
        "qc_fail_low_genes": adata.obs["n_genes"] < t.min_genes,
        "qc_fail_high_control": adata.obs["control_rate"] > t.max_control_rate,
    }
    if area is not None:
        flags["qc_fail_small_area"] = adata.obs["xa_cell_area"] < t.min_cell_area
        flags["qc_fail_large_area"] = adata.obs["xa_cell_area"] > t.max_cell_area
    if "nucleus_ratio" in adata.obs:
        flags["qc_fail_nucleus_ratio"] = adata.obs["nucleus_ratio"] > t.max_nucleus_ratio
    for name, values in flags.items():
        adata.obs[name] = np.asarray(values, dtype=bool)

    fail_cols = list(flags)
    adata.obs["qc_pass"] = ~adata.obs[fail_cols].to_numpy().any(axis=1)
    adata.uns["qc_thresholds"] = t.as_dict()
    adata.uns["qc_flag_columns"] = fail_cols
    return adata


def cell_qc_summary(adata: ad.AnnData) -> pd.Series:
    """One-row summary of the per-cell QC state of a sample."""
    obs = adata.obs
    out = {
        "sample": str(obs["sample"].iloc[0]) if "sample" in obs else "sample",
        "n_cells": int(adata.n_obs),
        "n_panel_genes": int(adata.n_vars),
        "median_transcripts_per_cell": float(obs["n_counts"].median()),
        "mean_transcripts_per_cell": float(obs["n_counts"].mean()),
        "median_genes_per_cell": float(obs["n_genes"].median()),
        "median_cell_area_um2": float(obs["xa_cell_area"].median()) if "xa_cell_area" in obs else np.nan,
        "median_nucleus_ratio": float(obs["nucleus_ratio"].median()) if "nucleus_ratio" in obs else np.nan,
        "frac_cells_no_nucleus": float(obs["no_nucleus"].mean()) if "no_nucleus" in obs else np.nan,
        "median_transcript_density": float(obs["transcript_density"].median()) if "transcript_density" in obs else np.nan,
        "frac_cells_passing": float(obs["qc_pass"].mean()),
    }
    for col in adata.uns.get("qc_flag_columns", []):
        out[f"frac_{col.replace('qc_fail_', 'fail_')}"] = float(obs[col].mean())
    return pd.Series(out)


def run_qc(adata: ad.AnnData, thresholds: QCThresholds | None = None) -> pd.Series:
    """Run-level background and yield metrics, with explicit pass/fail verdicts.

    The negative-control rates are computed from the data rather than copied
    from ``metrics_summary.csv`` so the number is reproducible from the matrix
    alone; the vendor value is carried alongside for comparison.
    """
    thresholds = thresholds or QCThresholds()
    obs = adata.obs
    panel_total = float(obs["n_counts"].sum())

    def _rate(col: str) -> float:
        return float(obs[col].sum()) / panel_total if col in obs and panel_total else 0.0

    neg_probe = _rate("negative_control_probe_counts")
    neg_codeword = _rate("negative_control_codeword_counts")
    unassigned = _rate("unassigned_codeword_counts")
    median_tx = float(obs["n_counts"].median())
    frac_pass = float(obs["qc_pass"].mean())

    checks = {
        "pass_neg_probe_rate": neg_probe <= thresholds.max_run_neg_probe_rate,
        "pass_neg_codeword_rate": neg_codeword <= thresholds.max_run_neg_codeword_rate,
        "pass_median_transcripts": median_tx >= thresholds.min_median_transcripts,
        "pass_cell_yield": frac_pass >= thresholds.min_cells_passing,
    }

    out = {
        "sample": str(obs["sample"].iloc[0]) if "sample" in obs else "sample",
        "n_cells": int(adata.n_obs),
        "total_panel_transcripts": int(panel_total),
        "neg_control_probe_rate": neg_probe,
        "neg_control_codeword_rate": neg_codeword,
        "unassigned_codeword_rate": unassigned,
        "median_transcripts_per_cell": median_tx,
        "frac_cells_passing": frac_pass,
        **{k: bool(v) for k, v in checks.items()},
        "run_verdict": "PASS" if all(checks.values()) else "REVIEW",
    }

    vendor = adata.uns.get("metrics_summary", {})
    for key in ("decoded_transcripts_per_100um2", "fraction_transcripts_assigned",
                "median_genes_per_cell", "median_transcripts_per_cell",
                "region_area", "num_cells_detected"):
        if key in vendor:
            out[f"vendor_{key}"] = vendor[key]
    return pd.Series(out)


# ---------------------------------------------------------------------------
# spatial QC
# ---------------------------------------------------------------------------

def grid_bin(adata: ad.AnnData, bin_um: float = 100.0) -> pd.DataFrame:
    """Aggregate per-cell QC onto a square grid of ``bin_um`` micrometres.

    Returns one row per occupied bin with cell density and the median of each
    QC metric, which is what makes regional artefacts visible.
    """
    xy = adata.obsm["spatial"]
    ix = np.floor((xy[:, 0] - xy[:, 0].min()) / bin_um).astype(int)
    iy = np.floor((xy[:, 1] - xy[:, 1].min()) / bin_um).astype(int)

    df = pd.DataFrame({
        "ix": ix, "iy": iy,
        "n_counts": adata.obs["n_counts"].to_numpy(),
        "control_rate": adata.obs["control_rate"].to_numpy(),
        "qc_pass": adata.obs["qc_pass"].to_numpy().astype(float),
    })
    if "xa_cell_area" in adata.obs:
        df["cell_area"] = adata.obs["xa_cell_area"].to_numpy()

    agg = df.groupby(["ix", "iy"]).agg(
        n_cells=("n_counts", "size"),
        median_counts=("n_counts", "median"),
        median_control_rate=("control_rate", "median"),
        frac_qc_pass=("qc_pass", "mean"),
        median_cell_area=("cell_area", "median") if "cell_area" in df else ("n_counts", "size"),
    ).reset_index()
    agg["x_um"] = agg["ix"] * bin_um + xy[:, 0].min()
    agg["y_um"] = agg["iy"] * bin_um + xy[:, 1].min()
    agg["cell_density_per_mm2"] = agg["n_cells"] / (bin_um ** 2) * 1e6
    return agg


def morans_i_grid(grid: pd.DataFrame, value_col: str, n_perm: int = 199,
                  seed: int = 0) -> dict:
    """Moran's I of a binned QC metric using rook contiguity on the grid.

    A metric that is spatially random across tissue has I near 0. Strong
    positive I means the metric comes in patches, which for a *technical*
    metric (control rate, cell density) is the signature of a regional
    artefact rather than of biology.
    """
    rng = np.random.default_rng(seed)
    coords = grid[["ix", "iy"]].to_numpy()
    values = grid[value_col].to_numpy(dtype=float)
    ok = np.isfinite(values)
    coords, values = coords[ok], values[ok]
    n = len(values)
    if n < 10:
        return {"metric": value_col, "morans_i": np.nan, "p_perm": np.nan, "n_bins": n}

    index = {(int(a), int(b)): i for i, (a, b) in enumerate(coords)}
    rows, cols = [], []
    for (a, b), i in index.items():
        for da, db in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            j = index.get((a + da, b + db))
            if j is not None:
                rows.append(i)
                cols.append(j)
    if not rows:
        return {"metric": value_col, "morans_i": np.nan, "p_perm": np.nan, "n_bins": n}
    W = sp.coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n)).tocsr()

    def _i(v: np.ndarray) -> float:
        z = v - v.mean()
        denom = (z ** 2).sum()
        if denom == 0:
            return np.nan
        return (n / W.sum()) * float(z @ (W @ z)) / denom

    observed = _i(values)
    null = np.array([_i(rng.permutation(values)) for _ in range(n_perm)])
    p = (np.sum(np.abs(null) >= abs(observed)) + 1) / (n_perm + 1)
    return {"metric": value_col, "morans_i": float(observed), "p_perm": float(p), "n_bins": int(n)}


def pseudobulk(adata: ad.AnnData, layer: str = "counts") -> pd.Series:
    """Summed panel-gene counts per sample, for cross-replicate comparison."""
    X = adata.layers[layer] if layer in adata.layers else adata.X
    totals = np.asarray(X.sum(axis=0)).ravel()
    return pd.Series(totals, index=adata.var_names, name=str(adata.obs["sample"].iloc[0]))
