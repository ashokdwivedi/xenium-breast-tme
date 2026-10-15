"""Readers for Xenium output bundles and the scRNA-seq reference atlas.

The Xenium ``outs`` bundle splits information that belongs together:

* ``cell_feature_matrix/``  - MEX triplet of the count matrix. Rows include the
  targeted panel genes *and* the negative-control features, which are the
  backbone of run-level QC and must not be dropped at read time.
* ``cells.csv.gz``          - one row per cell: centroid, segmented areas and
  the per-cell control counts.

:func:`read_xenium` joins the two into a single :class:`~anndata.AnnData` whose
``.X`` holds panel genes only, with every control feature preserved as a column
in ``.obs`` and the full feature table kept in ``.uns['control_features']``.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

# Feature classes emitted by the Xenium Analyzer. Everything that is not
# "Gene Expression" is a control used to estimate background, not biology.
PANEL_FEATURE_TYPE = "Gene Expression"
CONTROL_FEATURE_TYPES = (
    "Negative Control Probe",
    "Negative Control Codeword",
    "Unassigned Codeword",
    "Genomic Control",
    "Deprecated Codeword",
)


def _read_mex(path: Path) -> ad.AnnData:
    """Read a 10x MEX triplet (``matrix.mtx.gz``/``features``/``barcodes``)."""
    mtx = scipy.io.mmread(str(path / "matrix.mtx.gz")).T.tocsr().astype(np.float32)
    features = pd.read_csv(
        path / "features.tsv.gz", sep="\t", header=None,
        names=["gene_id", "gene_name", "feature_type"],
    )
    barcodes = pd.read_csv(path / "barcodes.tsv.gz", sep="\t", header=None, names=["cell_id"])

    adata = ad.AnnData(X=mtx)
    adata.var = features.set_index("gene_name")
    adata.var_names_make_unique()
    adata.obs_names = barcodes["cell_id"].astype(str).values
    return adata


def read_xenium(sample_dir: str | Path, sample_id: str | None = None) -> ad.AnnData:
    """Load one Xenium sample into an AnnData.

    Parameters
    ----------
    sample_dir
        Directory holding ``cell_feature_matrix/``, ``cells.csv.gz`` and
        (optionally) ``metrics_summary.csv`` and ``gene_panel.json``.
    sample_id
        Label written to ``.obs['sample']``. Defaults to the directory name.

    Returns
    -------
    AnnData
        ``.X``      raw integer counts for the panel genes only
        ``.obsm['spatial']``  (x, y) centroid in micrometres
        ``.obs``    segmentation areas plus one column per control feature class
        ``.uns``    ``metrics_summary``, ``panel``, ``control_features``
    """
    sample_dir = Path(sample_dir)
    sample_id = sample_id or sample_dir.name

    adata = _read_mex(sample_dir / "cell_feature_matrix")

    # --- split panel genes from control features -------------------------
    is_panel = adata.var["feature_type"].eq(PANEL_FEATURE_TYPE).values
    controls = adata[:, ~is_panel].copy()
    adata = adata[:, is_panel].copy()

    # Collapse each control class into a per-cell total.
    ctrl_counts = {}
    for ftype in controls.var["feature_type"].unique():
        cols = controls.var["feature_type"].eq(ftype).values
        key = ftype.lower().replace(" ", "_") + "_counts"
        ctrl_counts[key] = np.asarray(controls[:, cols].X.sum(axis=1)).ravel()
    for key, vals in ctrl_counts.items():
        adata.obs[key] = vals
    adata.uns["control_features"] = controls.var.reset_index()
    adata.uns["n_control_features"] = {
        str(k): int(v) for k, v in controls.var["feature_type"].value_counts().items()
    }

    # --- per-cell metadata -----------------------------------------------
    cells = pd.read_csv(sample_dir / "cells.csv.gz")
    cells["cell_id"] = cells["cell_id"].astype(str)
    cells = cells.set_index("cell_id").reindex(adata.obs_names)
    if cells["x_centroid"].isna().any():
        raise ValueError(f"{sample_id}: cells.csv.gz does not cover every matrix barcode")

    adata.obsm["spatial"] = cells[["x_centroid", "y_centroid"]].to_numpy(dtype=np.float64)
    for col in ("cell_area", "nucleus_area", "transcript_counts",
                "control_probe_counts", "control_codeword_counts", "total_counts"):
        if col in cells.columns:
            adata.obs[f"xa_{col}"] = cells[col].to_numpy()

    adata.obs["sample"] = pd.Categorical([sample_id] * adata.n_obs)
    adata.var["feature_type"] = adata.var["feature_type"].astype("category")

    # --- run-level context -------------------------------------------------
    metrics = sample_dir / "metrics_summary.csv"
    if metrics.exists():
        adata.uns["metrics_summary"] = pd.read_csv(metrics).iloc[0].to_dict()
    panel = sample_dir / "gene_panel.json"
    if panel.exists():
        with open(panel) as fh:
            meta = json.load(fh).get("metadata", {})
        adata.uns["panel"] = meta

    adata.layers["counts"] = adata.X.copy()
    return adata


def panel_gene_names(sample_dir: str | Path) -> list[str]:
    """Panel gene symbols for a Xenium sample, read from the features table."""
    features = pd.read_csv(
        Path(sample_dir) / "cell_feature_matrix" / "features.tsv.gz",
        sep="\t", header=None, names=["gene_id", "gene_name", "feature_type"],
    )
    return features.loc[features["feature_type"].eq(PANEL_FEATURE_TYPE), "gene_name"].tolist()


# ---------------------------------------------------------------------------
# scRNA-seq reference (Wu et al. 2021, GSE176078)
# ---------------------------------------------------------------------------

def read_reference_subset(
    ref_dir: str | Path,
    keep_genes: list[str],
    chunk_rows: int = 5_000_000,
) -> ad.AnnData:
    """Stream the Wu et al. atlas, keeping only ``keep_genes``.

    The full matrix is ~100k cells x ~29k genes and does not need to be
    materialised: label transfer only ever uses the Xenium panel genes. The
    MatrixMarket body is parsed in chunks and filtered on the fly, which keeps
    peak memory at a few hundred MB instead of several GB.
    """
    ref_dir = Path(ref_dir)
    genes = pd.read_csv(ref_dir / "count_matrix_genes.tsv", header=None)[0].astype(str).tolist()
    barcodes = pd.read_csv(ref_dir / "count_matrix_barcodes.tsv", header=None)[0].astype(str).tolist()

    wanted = {g: i for i, g in enumerate(keep_genes)}
    # 1-based MatrixMarket row index -> position in the output matrix
    row_map = np.full(len(genes) + 1, -1, dtype=np.int32)
    for i, g in enumerate(genes):
        if g in wanted:
            row_map[i + 1] = wanted[g]

    mtx_path = ref_dir / "count_matrix_sparse.mtx"
    opener = gzip.open if mtx_path.suffix == ".gz" else open
    with opener(mtx_path, "rt") as fh:  # type: ignore[operator]
        header = fh.readline()
        if not header.startswith("%%MatrixMarket"):
            raise ValueError("not a MatrixMarket file")
        line = fh.readline()
        while line.startswith("%"):
            line = fh.readline()
        n_genes, n_cells, _nnz = (int(v) for v in line.split())
        if n_genes != len(genes) or n_cells != len(barcodes):
            raise ValueError("mtx dimensions disagree with the genes/barcodes files")

        rows, cols, vals = [], [], []
        reader = pd.read_csv(
            fh, sep=r"\s+", header=None, names=["g", "c", "v"],
            dtype={"g": np.int32, "c": np.int32, "v": np.float32}, chunksize=chunk_rows,
        )
        for chunk in reader:
            mapped = row_map[chunk["g"].to_numpy()]
            keep = mapped >= 0
            if keep.any():
                rows.append(mapped[keep])
                cols.append(chunk["c"].to_numpy()[keep] - 1)
                vals.append(chunk["v"].to_numpy()[keep])

    X = sp.coo_matrix(
        (np.concatenate(vals), (np.concatenate(cols), np.concatenate(rows))),
        shape=(n_cells, len(keep_genes)),
    ).tocsr().astype(np.float32)

    adata = ad.AnnData(X=X)
    adata.obs_names = barcodes
    adata.var_names = keep_genes

    meta = pd.read_csv(ref_dir / "metadata.csv", index_col=0)
    meta.index = meta.index.astype(str)
    adata.obs = adata.obs.join(meta.reindex(adata.obs_names))
    adata.layers["counts"] = adata.X.copy()
    return adata
