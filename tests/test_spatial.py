import numpy as np
import pytest

from xtme import annotate, spatial


def test_neighbourhood_composition_rows_sum_to_one(toy_typed):
    comp = spatial.neighbourhood_composition(toy_typed, "cell_type", k=10)
    np.testing.assert_allclose(comp.sum(axis=1).to_numpy(), 1.0, atol=1e-6)
    assert list(comp.columns) == ["TypeA", "TypeB"]
    assert len(comp) == toy_typed.n_obs


def test_neighbourhood_composition_recovers_compartments(toy_typed):
    """Cells in a spatially pure compartment see mostly their own type."""
    comp = spatial.neighbourhood_composition(toy_typed, "cell_type", k=10)
    assert comp.iloc[:300]["TypeA"].mean() > 0.9
    assert comp.iloc[300:]["TypeB"].mean() > 0.9


def test_niches_separate_the_two_compartments(toy_typed):
    niche, profile = spatial.compute_niches(toy_typed, "cell_type", k=10, n_niches=2, seed=0)
    assert len(niche) == toy_typed.n_obs
    assert profile.shape == (2, 2)
    np.testing.assert_allclose(profile.sum(axis=1).to_numpy(), 1.0, atol=1e-6)
    # each niche should be dominated by a different cell type
    assert set(profile.idxmax(axis=1)) == {"TypeA", "TypeB"}


def test_name_niches_mentions_top_type(toy_typed):
    _, profile = spatial.compute_niches(toy_typed, "cell_type", k=10, n_niches=2, seed=0)
    names = spatial.name_niches(profile, top_n=1)
    assert len(names) == 2
    assert all(("TypeA" in v) or ("TypeB" in v) for v in names.values())


def test_distance_to_compartment_is_zero_inside(toy_typed):
    mask = toy_typed.obs["cell_type"].eq("TypeA").to_numpy()
    d = spatial.distance_to_compartment(toy_typed, mask)
    assert d.shape == (toy_typed.n_obs,)
    np.testing.assert_allclose(d[mask], 0.0, atol=1e-9)
    # the other compartment is physically separated by >=200 um by construction
    assert d[~mask].min() > 150


def test_distance_to_compartment_rejects_empty_mask(toy_typed):
    with pytest.raises(ValueError):
        spatial.distance_to_compartment(toy_typed, np.zeros(toy_typed.n_obs, dtype=bool))


def test_interface_enrichment_null_is_centred(toy_typed):
    """With a single non-target class the observed value must equal the null."""
    toy_typed.obs["dist_to_tumour"] = spatial.distance_to_compartment(
        toy_typed, toy_typed.obs["cell_type"].eq("TypeA").to_numpy()
    )
    out = spatial.interface_enrichment(
        toy_typed, "cell_type", "TypeA", interface_um=1e6, n_perms=99, seed=0
    )
    assert len(out) == 1
    assert out.iloc[0]["frac_within_interface"] == pytest.approx(1.0)
    assert out.iloc[0]["log2_enrichment"] == pytest.approx(0.0, abs=1e-6)


def test_interface_enrichment_columns(toy_typed):
    toy_typed.obs["dist_to_tumour"] = spatial.distance_to_compartment(
        toy_typed, toy_typed.obs["cell_type"].eq("TypeA").to_numpy()
    )
    out = spatial.interface_enrichment(toy_typed, "cell_type", "TypeA",
                                       interface_um=300, n_perms=99, seed=0)
    for col in ("cell_type", "frac_within_interface", "null_mean", "log2_enrichment",
                "z_score", "p_perm", "q_bh"):
        assert col in out.columns
    assert ((out["p_perm"] > 0) & (out["p_perm"] <= 1)).all()
    assert (out["q_bh"] <= 1).all()


def test_distance_profile_rows_sum_to_one(toy_typed):
    toy_typed.obs["dist_to_tumour"] = spatial.distance_to_compartment(
        toy_typed, toy_typed.obs["cell_type"].eq("TypeA").to_numpy()
    )
    prof = spatial.distance_profile(toy_typed, "cell_type", "dist_to_tumour")
    occupied = prof.loc[prof.sum(axis=1) > 0]
    np.testing.assert_allclose(occupied.sum(axis=1).to_numpy(), 1.0, atol=1e-6)


def test_marker_labelling_flags_ambiguity():
    import pandas as pd
    scores = pd.DataFrame(
        {"T cell": [2.0, 0.50], "Myeloid": [0.1, 0.49]}, index=["0", "1"]
    )
    out = annotate.label_clusters_by_marker(scores, margin=0.10)
    assert out.loc["0", "marker_label"] == "T cell"
    assert not out.loc["0", "ambiguous"]
    assert out.loc["1", "ambiguous"]


def test_concordance_rows_normalise():
    import pandas as pd
    a = pd.Series(["x", "x", "y", "y"])
    b = pd.Series(["p", "q", "q", "q"])
    tab = annotate.concordance(a, b)
    np.testing.assert_allclose(tab.sum(axis=1).to_numpy(), 1.0)
    assert tab.loc["y", "q"] == pytest.approx(1.0)
