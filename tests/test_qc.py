import numpy as np
import pytest

from xtme import qc


def test_cell_qc_adds_expected_columns(toy_xenium):
    qc.compute_cell_qc(toy_xenium)
    for col in ("n_counts", "n_genes", "control_counts", "control_rate",
                "nucleus_ratio", "transcript_density", "qc_pass"):
        assert col in toy_xenium.obs


def test_counts_match_the_matrix(toy_xenium):
    qc.compute_cell_qc(toy_xenium)
    expected = np.asarray(toy_xenium.layers["counts"].sum(axis=1)).ravel()
    np.testing.assert_allclose(toy_xenium.obs["n_counts"].to_numpy(), expected)


def test_degenerate_cells_are_flagged(toy_xenium):
    """The 10 seeded bad cells must fail, and they must be the only area failures."""
    qc.compute_cell_qc(toy_xenium)
    obs = toy_xenium.obs
    bad = obs.index[:10]
    assert obs.loc[bad, "qc_pass"].sum() == 0
    assert obs.loc[bad, "qc_fail_low_counts"].all()
    assert obs.loc[bad, "qc_fail_small_area"].all()
    assert obs["qc_fail_small_area"].sum() == 10


def test_thresholds_are_respected(toy_xenium):
    strict = qc.QCThresholds(min_transcripts=10_000)
    qc.compute_cell_qc(toy_xenium, strict)
    assert toy_xenium.obs["qc_pass"].sum() == 0
    assert toy_xenium.uns["qc_thresholds"]["min_transcripts"] == 10_000


def test_control_rate_bounds(toy_xenium):
    qc.compute_cell_qc(toy_xenium)
    rate = toy_xenium.obs["control_rate"].to_numpy()
    assert np.isfinite(rate).all()
    assert ((rate >= 0) & (rate <= 1)).all()


def test_run_qc_verdict_and_rates(toy_xenium):
    qc.compute_cell_qc(toy_xenium)
    row = qc.run_qc(toy_xenium)
    assert row["run_verdict"] in {"PASS", "REVIEW"}
    assert 0 <= row["neg_control_probe_rate"] <= 1
    assert row["n_cells"] == toy_xenium.n_obs
    # the verdict must be the conjunction of its own checks
    checks = [row[k] for k in row.index if k.startswith("pass_")]
    assert (row["run_verdict"] == "PASS") == all(checks)


def test_grid_bin_covers_every_cell(toy_xenium):
    qc.compute_cell_qc(toy_xenium)
    grid = qc.grid_bin(toy_xenium, bin_um=100.0)
    assert grid["n_cells"].sum() == toy_xenium.n_obs
    assert (grid["cell_density_per_mm2"] > 0).all()


def test_morans_i_detects_structure(toy_xenium):
    """A metric painted onto one half of the section must be autocorrelated."""
    qc.compute_cell_qc(toy_xenium)
    grid = qc.grid_bin(toy_xenium, bin_um=100.0)
    grid["left_half"] = (grid["x_um"] < grid["x_um"].median()).astype(float)
    structured = qc.morans_i_grid(grid, "left_half", n_perm=199, seed=0)
    assert structured["morans_i"] > 0.3
    assert structured["p_perm"] < 0.05


def test_morans_i_on_noise_is_near_zero(toy_xenium):
    rng = np.random.default_rng(1)
    qc.compute_cell_qc(toy_xenium)
    grid = qc.grid_bin(toy_xenium, bin_um=100.0)
    grid["noise"] = rng.normal(size=len(grid))
    out = qc.morans_i_grid(grid, "noise", n_perm=199, seed=0)
    assert abs(out["morans_i"]) < 0.3


def test_pseudobulk_totals(toy_xenium):
    pb = qc.pseudobulk(toy_xenium)
    assert pb.sum() == pytest.approx(toy_xenium.layers["counts"].sum())
    assert list(pb.index) == list(toy_xenium.var_names)
