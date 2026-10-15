#!/usr/bin/env python
"""Stage 4 - build the project website into ``docs/``.

Every figure is copied from ``results/figures`` and every number in the prose
is read from ``results/tables`` at build time, so the site cannot drift away
from the analysis that produced it.

    python scripts/04_build_site.py [--repo-url URL]
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from sitelib import (  # noqa: E402
    AUTHOR_SITE,
    PAGES,
    cards,
    caveat,
    figure,
    note,
    stats,
    table,
    write,
)

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from xtme.config import load_config  # noqa: E402

CFG = load_config()

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "results" / "figures"
TAB = ROOT / "results" / "tables"
DOCS = ROOT / "docs"
DEFAULT_REPO = "https://github.com/ashokdwivedi/xenium-breast-tme"


def fmt(x: float, nd: int = 2) -> str:
    return f"{x:,.{nd}f}"


def pct(x: float, nd: int = 1) -> str:
    return f"{x * 100:.{nd}f}%"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-url", default=DEFAULT_REPO)
    args = ap.parse_args()
    repo = args.repo_url.rstrip("/")

    DOCS.mkdir(exist_ok=True)
    figdir = DOCS / "assets" / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    for png in sorted(FIG.glob("*.png")):
        shutil.copy2(png, figdir / png.name)

    # also ship the tables so every claim is checkable without a clone
    tabdir = DOCS / "assets" / "tables"
    tabdir.mkdir(parents=True, exist_ok=True)
    for csv in sorted(TAB.glob("*.csv")):
        if csv.stat().st_size < 2_000_000:
            shutil.copy2(csv, tabdir / csv.name)

    d = load_results()
    write(DOCS, PAGES[0], page_overview(d, repo), repo)
    write(DOCS, PAGES[1], page_qc(d), repo)
    write(DOCS, PAGES[2], page_celltypes(d), repo)
    write(DOCS, PAGES[3], page_spatial(d), repo)
    write(DOCS, PAGES[4], page_methods(d, repo), repo)
    write(DOCS, PAGES[5], page_about(repo), repo)
    (DOCS / ".nojekyll").write_text("")
    print(f"[04] wrote {len(PAGES)} pages to {DOCS}")


# ---------------------------------------------------------------------------
# results loading
# ---------------------------------------------------------------------------
def load_results() -> dict:
    sc = pd.read_csv(TAB / "qc_scorecard.csv", index_col=0)
    thr = pd.read_csv(TAB / "qc_thresholds.csv", index_col=0)["value"]
    moran = pd.read_csv(TAB / "spatial_autocorrelation.csv")
    rep = pd.read_csv(TAB / "replicate_summary.csv", index_col=0)["value"]
    conc = pd.read_csv(TAB / "replicate_gene_concordance.csv", index_col=0)
    cts = pd.read_csv(TAB / "celltype_counts.csv", index_col=0)
    ann = json.loads((TAB / "annotation_summary.json").read_text())
    clus = pd.read_csv(TAB / "cluster_annotation_table.csv", index_col=0)
    cov = pd.read_csv(TAB / "marker_coverage.csv")
    ne = pd.read_csv(TAB / "neighborhood_enrichment.csv", index_col=0)
    pairs = pd.read_csv(TAB / "neighborhood_pairs.csv")
    niche = pd.read_csv(TAB / "niche_composition.csv", index_col=0)
    iface = pd.read_csv(TAB / "immune_distance_stats.csv")
    return dict(sc=sc, thr=thr, moran=moran, rep=rep, conc=conc, cts=cts, ann=ann,
                clus=clus, cov=cov, ne=ne, pairs=pairs, niche=niche, iface=iface)


# ---------------------------------------------------------------------------
# pages
# ---------------------------------------------------------------------------
def page_overview(d: dict, repo: str) -> str:
    sc, ann, iface = d["sc"], d["ann"], d["iface"]
    total_cells = int(sc["n_cells"].sum())
    enr = iface[iface["cell_type"] != "Cancer Epithelial"].sort_values(
        "log2_enrichment", ascending=False)
    top_in = enr.iloc[0]
    top_out = enr.iloc[-1]

    body = f"""
<p>This is a complete, reproducible analysis of a public <b>10x Genomics Xenium In Situ</b>
run on FFPE human breast cancer: two replicate sections of the same tumour block, taken from
raw vendor output through quality control, cell typing against a single-cell reference, and
spatial statistics on how the tumour and immune compartments are arranged relative to each
other. Everything here \u2014 the code, the thresholds, the unit tests, the figures and the
tables behind every number on this site \u2014 is in
<a href="{repo}">the repository</a>, and the whole thing reruns from a clean clone with
<code>bash scripts/download_data.sh &amp;&amp; make all</code>.</p>

{stats([
    (f"{total_cells:,}", "segmented cells across two replicate sections"),
    (f"{int(sc['n_panel_genes'].iloc[0])}", "panel genes (targeted, not whole-transcriptome)"),
    (pct(sc['frac_cells_passing'].min(), 1), "cells passing QC in the stricter replicate"),
    (fmt(d['rep']['pearson_log10cpm'], 3), "cross-replicate pseudobulk correlation (r)"),
    (f"{ann['n_leiden_clusters']}", "transcriptional clusters typed against the reference"),
    (f"{len(d['niche'])}", "recurrent spatial niches recovered"),
])}

<h2>What the analysis shows</h2>

<p>Three findings carry the story, and each is a link to the page that derives it.</p>

<ol>
<li><b>Both runs are technically sound, and the evidence is quantitative rather than
    visual.</b> Negative-control probe rates are
    {fmt(sc['neg_control_probe_rate'].min() * 100, 3)}\u2013{fmt(sc['neg_control_probe_rate'].max() * 100, 3)}%
    of panel signal, median depth is identical at
    {int(sc['median_transcripts_per_cell'].iloc[0])} transcripts per cell in both sections, and
    pseudobulk expression across the {int(d['rep']['n_genes'])} panel genes correlates at
    r&nbsp;=&nbsp;{fmt(d['rep']['pearson_log10cpm'], 3)} between replicates. The spatial
    structure that <em>is</em> present in the QC metrics tracks tissue architecture, not
    instrument geometry \u2014 which is the distinction that decides whether you filter cells
    or re-run the slide. <a href="qc.html">See the QC page.</a></li>

<li><b>Cell typing on a {int(sc['n_panel_genes'].iloc[0])}-gene panel is reliable for broad
    lineages and honest about where it stops.</b> Label transfer from the Wu&nbsp;et&nbsp;al.
    breast-cancer atlas assigns {pct(1 - ann['frac_unassigned'], 1)} of QC-passing cells above
    a confidence floor, and an independent marker-scoring route agrees with it on
    {pct(ann['frac_clusters_routes_agree'], 0)} of clusters. The disagreements are reported,
    not smoothed over. <a href="celltypes.html">See the cell-typing page.</a></li>

<li><b>The tissue is compartmentalised, and the tumour boundary is selectively permeable.</b>
    Against a label-permutation null that holds every cell position fixed,
    <b>{top_in['cell_type']}</b> is the population most enriched within
    {int(CFG['spatial']['interface_um'])}\u00b5m of the tumour compartment
    (log<sub>2</sub> {fmt(top_in['log2_enrichment'])}), while <b>{top_out['cell_type']}</b> is
    the most excluded (log<sub>2</sub> {fmt(top_out['log2_enrichment'])}).
    <a href="spatial.html">See the spatial page.</a></li>
</ol>

<h2>Read it in order</h2>

{cards([
    ("Stage 1", "qc.html", "Quality control",
     "Per-cell, per-run and position-aware checks, with the pass/fail criteria stated up front "
     "and a replicate-concordance arm that separates biology from batch."),
    ("Stage 2", "celltypes.html", "Cell typing",
     "Harmony-integrated k-NN label transfer from an scRNA-seq atlas, cross-checked against "
     "canonical markers, with panel coverage reported as a limit on resolution."),
    ("Stage 3", "spatial.html", "Spatial structure",
     "Neighbourhood enrichment, recurrent niches, and permutation-tested enrichment at the "
     "tumour\u2013immune interface."),
    ("Reference", "methods.html", "Methods and limits",
     "Every parameter and threshold, the reproduction recipe, and a frank list of what this "
     "analysis cannot establish."),
])}

<h2>Why it was built this way</h2>

<p>Three choices distinguish this from a notebook walkthrough, and each one exists because
the alternative fails in production.</p>

<p><b>Thresholds are configuration, not folklore.</b> Every QC cutoff lives in
<code>config/config.yaml</code> and is written into the output tables alongside the verdict it
produced. When a criterion is later judged too strict, the change is a one-line diff with a
visible blast radius, not an archaeology exercise across notebook cells.</p>

<p><b>Both replicates are analysed, not just the prettier one.</b> A single section cannot
distinguish a biological gradient from a slide-level artefact. Two sections of the same block
can, and the cross-replicate comparison is reported as a first-class result rather than a
supplementary figure.</p>

<p><b>Cell typing is run twice, by independent routes, and the disagreement is published.</b>
Reference transfer and marker scoring fail in different ways: transfer inherits the reference's
composition bias, marker scoring inherits the panel's blind spots. Where they agree the call is
solid; where they do not, the honest output is a flag, not a label.</p>

{note(f'''<p>The repository is tested. {_test_count()} unit tests run on synthetic fixtures in
CI on every push \u2014 no data download required \u2014 covering the QC flags, the Moran's I
statistic, niche recovery, the interface null model and the marker-ambiguity rule. A
statistical function that cannot be tested on data with a known answer is a function nobody
should trust.</p>''')}

<p style="color:var(--ink-faint);font-size:0.9rem">Analysis and site by
<a href="{AUTHOR_SITE}">Ashok Dwivedi</a> \u2014 <a href="{AUTHOR_SITE}">{AUTHOR_SITE.replace('https://', '')}</a>.</p>
"""
    return body


def _test_count() -> str:
    """Count test functions so the claim on the site stays true."""
    n = sum(
        line.strip().startswith("def test_")
        for f in (ROOT / "tests").glob("test_*.py")
        for line in f.read_text().splitlines()
    )
    return str(n)


def page_qc(d: dict) -> str:
    sc, thr, moran = d["sc"], d["thr"], d["moran"]
    r1, r2 = sc.index[0], sc.index[1]
    rows = []
    label_map = [
        ("n_cells", "Segmented cells", 0, "{:,.0f}"),
        ("median_transcripts_per_cell", "Median transcripts / cell", 0, "{:,.0f}"),
        ("median_genes_per_cell", "Median genes / cell", 0, "{:,.0f}"),
        ("median_cell_area_um2", "Median cell area (\u00b5m\u00b2)", 1, "{:,.1f}"),
        ("neg_control_probe_rate", "Negative-control probe rate", 5, "{:.5f}"),
        ("neg_control_codeword_rate", "Negative-control codeword rate", 6, "{:.6f}"),
        ("vendor_fraction_transcripts_assigned", "Transcripts assigned to a cell (vendor)", 3, "{:.3f}"),
        ("frac_cells_passing", "Cells passing all per-cell criteria", 4, "{:.4f}"),
    ]
    for key, label, _nd, f in label_map:
        rows.append([label, f.format(sc.loc[r1, key]), f.format(sc.loc[r2, key])])
    verdict = [f'<span class="pass">{sc.loc[s, "run_verdict"]}</span>' for s in (r1, r2)]
    rows.append(["<b>Run verdict</b>", f"<b>{verdict[0]}</b>", f"<b>{verdict[1]}</b>"])

    fail_rows = [
        [lbl, pct(sc.loc[r1, k], 2), pct(sc.loc[r2, k], 2)]
        for k, lbl in [
            ("frac_fail_low_counts", f"Fewer than {int(thr['min_transcripts'])} transcripts"),
            ("frac_fail_low_genes", f"Fewer than {int(thr['min_genes'])} genes detected"),
            ("frac_fail_small_area", f"Cell area &lt; {thr['min_cell_area']:.0f} \u00b5m\u00b2"),
            ("frac_fail_large_area", f"Cell area &gt; {thr['max_cell_area']:.0f} \u00b5m\u00b2"),
            ("frac_fail_high_control", f"Control fraction &gt; {thr['max_control_rate']:.0%}"),
            ("frac_fail_nucleus_ratio", "Nucleus larger than its cell"),
        ]
    ]

    mo = moran.pivot(index="metric", columns="sample", values="morans_i")
    mo_rows = [[_pretty_metric(m), fmt(mo.loc[m, r1], 3), fmt(mo.loc[m, r2], 3)]
               for m in ["cell_density_per_mm2", "median_counts", "frac_qc_pass",
                         "median_control_rate"]]
    dens_i = float(mo.loc["cell_density_per_mm2"].max())
    ctrl_i = float(mo.loc["median_control_rate"].min())

    body = f"""
<p>Quality control for in situ data is not one question but three, and they are worth keeping
apart because they fail for different reasons and have different remedies. <b>Is this cell
real?</b> is a filtering problem. <b>Is this run usable?</b> is a go/no-go decision about the
slide. <b>Is the defect regional?</b> is the question that distinguishes a tissue gradient from
a decoding or imaging artefact, and it is the one most often skipped.</p>

<h2>The criteria, stated before the data</h2>

<p>Every threshold below is read from <code>config/config.yaml</code>, written into
<code>results/tables/qc_thresholds.csv</code> next to the verdict it produced, and applied
identically to both sections. Nothing was tuned after seeing the outcome.</p>

{table(["Per-cell criterion", "Cells failing, " + str(r1), "Cells failing, " + str(r2)], fail_rows)}

<p>The criteria are deliberately permissive. On a targeted panel the dominant failure mode is
not a noisy cell but an <em>empty</em> one \u2014 a segmentation mask with no transcripts in it
\u2014 and a conservative depth floor removes those without discarding genuinely small immune
cells, which is the error that quietly deletes the most interesting population in a tumour
section.</p>

<h2>Per-cell distributions</h2>

{figure("assets/figures/fig01_qc_panel.png",
        "<b>Per-cell quality across both replicates.</b> (a) Depth is near-identical between "
        "sections over four orders of magnitude; the knee at the right edge is the empty-mask "
        "population the depth filter is designed to remove. (b, c) Transcript and gene "
        "distributions are unimodal above the filter, with the spike at one transcript per cell "
        "clearly separated from the main mode. (d) Control-feature background concentrates "
        "almost entirely in shallow cells \u2014 the diagonal band is the arithmetic floor of a "
        "single control count, not a gradient of contamination. (e) Nucleus-to-cell area ratio "
        "stays below one everywhere, so segmentation expansion never produced an impossible "
        "cell. (f) Which criterion actually removes cells, ranked.",
        alt="Six-panel figure of per-cell QC distributions for both replicates")}

<h2>Run-level summary</h2>

<p>Run metrics are recomputed from the count matrix rather than copied from the vendor summary,
with the vendor value carried alongside so the two can be compared. The negative-control rate
is the single most informative number on this table: it is the fraction of detected signal that
is known to be false, and it bounds how much of any weak biological effect could be background.</p>

{table(["Metric", str(r1), str(r2)], rows)}

{note(f'''<p>Both runs pass every criterion. Negative-control probes account for
{pct(sc.loc[r1, "neg_control_probe_rate"], 3)} and {pct(sc.loc[r2, "neg_control_probe_rate"], 3)}
of panel signal respectively \u2014 roughly one false transcript in a thousand \u2014 and
negative-control <em>codewords</em>, which report decoding error rather than probe
non-specificity, are an order of magnitude rarer still. The cell yield differs substantially
between sections ({sc.loc[r1, "n_cells"]:,.0f} versus {sc.loc[r2, "n_cells"]:,.0f}), but this
tracks the amount of tissue on the slide rather than any difference in run quality: median depth
and median panel breadth per cell are identical.</p>''')}

<h2>Is any defect regional?</h2>

<p>A run-level average can hide a bad corner. The test here is spatial autocorrelation: bin the
section into {int(d['thr']['grid_bin_um'])}&nbsp;\u00b5m squares, compute Moran's I for each
technical metric on the queen-contiguity grid, and ask whether the metric clusters in space more
than a random relabelling of the same bins would produce. A technical artefact \u2014 a defocused
tile, a decoding dropout, a fluidics streak \u2014 produces strong, blocky autocorrelation aligned
to instrument geometry. Tissue produces strong autocorrelation too, but shaped like tissue.</p>

{table(["Metric", f"Moran's I, {r1}", f"Moran's I, {r2}"], mo_rows)}

{figure("assets/figures/fig02_spatial_qc.png",
        "<b>Technical metrics mapped back onto the section.</b> Depth (left) and cell density "
        "(right) are strongly structured, but the structure is ductal and stromal architecture: "
        "the dense, transcript-rich regions are epithelial nests with sharp, curved boundaries, "
        "not rectangular fields of view. Control-feature background (centre) is close to flat. "
        "Colour scales are shared between rows so the two sections are directly comparable.",
        alt="Spatial maps of transcript depth, control background and cell density for both replicates")}

<p>Cell density is the most strongly autocorrelated metric (I&nbsp;=&nbsp;{fmt(dens_i, 2)}),
which is exactly what a tissue with dense epithelial nests and sparse stroma should produce.
Control-feature background is the least structured (I&nbsp;=&nbsp;{fmt(ctrl_i, 2)} in
{r2}) \u2014 the signature of background that is uniform across the slide. Crucially, none of
the maps show rectangular blocks or straight-edged discontinuities. The depth structure follows
curved tissue boundaries, so the correct action is to account for cell type when interpreting
depth, not to mask a region or re-run the slide.</p>

{caveat('''<p>Moran's I on a binned grid detects regional structure; it does not attribute it.
A defect perfectly confounded with a tissue boundary would be indistinguishable from biology by
this test alone. Registering the DAPI image against an H&amp;E section and inspecting the
flagged bins directly is the next check, and it needs the morphology image that is not part of
this public bundle.</p>''')}

<h2>Do the two replicates agree?</h2>

<p>This is the arm that a single-section analysis cannot have. Two serial sections of the same
block were processed as independent runs, so any difference between them is technical by
construction \u2014 sectioning, hybridisation, imaging, decoding \u2014 and sets a floor on the
effect size worth believing from one section.</p>

{figure("assets/figures/fig03_replicate_concordance.png",
        f"<b>Cross-replicate concordance of panel-level expression.</b> (a) Pseudobulk CPM for "
        f"all {int(d['rep']['n_genes'])} panel genes, r&nbsp;=&nbsp;{fmt(d['rep']['pearson_log10cpm'], 3)} "
        f"on log\u2081\u2080 values. (b) The same data as an MA plot: only "
        f"{int(d['rep']['n_genes_abs_log2fc_gt1'])} genes exceed a two-fold difference, and they "
        f"are scattered across the abundance range rather than concentrated at low expression, "
        f"which is what a depth or sensitivity difference would produce.",
        alt="Scatter and MA plots comparing pseudobulk expression between the two replicates")}

<p>Agreement is high (r&nbsp;=&nbsp;{fmt(d['rep']['pearson_log10cpm'], 3)},
\u03c1&nbsp;=&nbsp;{fmt(d['rep']['spearman_cpm'], 3)}), and the
{int(d['rep']['n_genes_abs_log2fc_gt1'])} genes that do differ by more than two-fold are
informative: the largest discordances are in {_top_genes(d['conc'])}. These are markers of
adipose and skeletal-muscle tissue \u2014 structures that are present in one section plane and
absent from the other. That is a sectioning difference, not a run-quality difference, and it is
the clearest possible demonstration of why a fold-change from a single section should not be
read as biology without a replicate to anchor it.</p>

{note('''<p><b>The practical output of this stage</b> is a number to calibrate against: on this
tissue, a two-fold pseudobulk difference between sections of the <em>same block</em> is
achievable from sectioning alone. Any cross-sample claim smaller than that needs more than one
section per sample to support it.</p>''')}
"""
    return body


def _pretty_metric(m: str) -> str:
    return {
        "cell_density_per_mm2": "Cell density (cells/mm\u00b2)",
        "median_counts": "Median transcripts per cell",
        "frac_qc_pass": "Fraction of cells passing QC",
        "median_control_rate": "Median control-feature fraction",
    }.get(m, m)


def _top_genes(conc: pd.DataFrame, n: int = 4) -> str:
    g = list(conc.sort_values("abs_log2fc", ascending=False).head(n).index)
    return ", ".join(f"<i>{x}</i>" for x in g[:-1]) + f" and <i>{g[-1]}</i>"


def page_celltypes(d: dict) -> str:
    cts, ann, clus, cov = d["cts"], d["ann"], d["clus"], d["cov"]
    rows = [[idx, f"{int(r['n_cells']):,}", pct(r["pct"] / 100, 1)]
            for idx, r in cts.iterrows()]
    missing = cov[cov["n_found"] < cov["n_requested"]]
    disagree = clus[~clus["routes_agree"]]

    body = f"""
<p>A {int(d['sc']['n_panel_genes'].iloc[0])}-gene targeted panel is not a transcriptome. It
resolves the lineages it was designed to resolve and is blind to everything else, so the useful
question is not "what cell type is this?" but "which assignments does this panel actually
support, and which are guesses dressed as labels?" The approach here is to type the cells twice,
by methods that fail differently, and report where they disagree.</p>

<h2>Route A \u2014 transfer from a single-cell reference</h2>

<p>The reference is the breast-cancer atlas of Wu&nbsp;et&nbsp;al. (2021, GEO GSE176078):
{ann['n_reference_cells']:,} cells after subsetting to panel genes and capping each label at
{int(CFG['reference']['max_cells_per_label']):,} cells, so that an abundant population cannot win the
k-nearest-neighbour vote on prior alone. Reference and Xenium cells are normalised identically,
projected into a shared PCA space, and integrated with Harmony \u2014 because the dominant axis
of variation between an scRNA-seq dataset and an in situ panel is the assay, not the biology,
and a k-NN vote in unintegrated space would simply recover the batch.</p>

<p>Each Xenium cell takes a distance-weighted vote of its
{int(CFG['transfer']['k'])} nearest reference neighbours. A cell whose vote does not reach the
confidence floor is labelled <code>Unassigned</code> rather than forced into the nearest class
\u2014 {pct(ann['frac_unassigned'], 2)} of cells end up there, and keeping them visible is the
point: a transfer that assigns everything is not more accurate, it is less honest.</p>

<h2>Route B \u2014 canonical markers, scored on clusters</h2>

<p>Independently, cells are clustered with Leiden
({ann['n_leiden_clusters']} clusters at the chosen resolution) and each cluster is scored
against curated marker sets for the lineages a breast panel can see. Scoring is done on
clusters rather than single cells because a {int(d['sc']['n_panel_genes'].iloc[0])}-gene panel
at single-cell depth is too sparse for per-cell marker calls to be stable. A cluster whose best
and second-best lineage scores are within a narrow margin is flagged <code>ambiguous</code>
instead of labelled.</p>

{figure("assets/figures/fig04_umap_clusters.png",
        "<b>Unsupervised structure, before any labels.</b> (a) Leiden clusters in UMAP space. "
        "(b) The same clusters mapped back to tissue coordinates \u2014 the clusters are "
        "spatially coherent, which is a meaningful check: transcriptional clusters derived "
        "without any positional information should not be spatially random in a structured "
        "tissue, and clusters that <em>are</em> scattered uniformly are usually technical.",
        alt="UMAP of Leiden clusters beside the same clusters plotted in tissue coordinates")}

<h2>The resulting map</h2>

{figure("assets/figures/fig05_celltypes.png",
        "<b>Cell types in transcriptional and physical space.</b> Colour is bound to cell-type "
        "identity once and reused in every figure on this site. The tissue map (centre) shows "
        "the architecture the composition bar chart (right) cannot: epithelial tumour in "
        "discrete nests, stromal and immune populations filling the space between them.",
        alt="UMAP, tissue map and composition bar chart of annotated cell types")}

{table(["Cell type", "Cells", "Share"], rows)}

<h2>Where the two routes disagree</h2>

<p>The routes agree on {pct(ann['frac_clusters_routes_agree'], 0)} of clusters. Agreement
between a reference-based and a marker-based call is not proof \u2014 both can be wrong the same
way if the panel lacks a discriminating gene \u2014 but disagreement is a reliable warning.</p>

{figure("assets/figures/fig06_annotation_quality.png",
        "<b>How far to trust each label.</b> (a) Transfer confidence by assigned type; the "
        "dashed line is the floor below which a cell is left unassigned. Populations with "
        "systematically low confidence are the ones the panel separates poorly. (b) Agreement "
        "matrix between the two routes, row-normalised. (c) Marker coverage per lineage \u2014 "
        "dark bars are canonical markers present on the panel, light bars are markers the panel "
        "does not carry. A lineage with few panel markers cannot be resolved by either route, "
        "and reading this panel first prevents over-interpreting the other two.",
        alt="Transfer confidence, route agreement matrix and marker coverage per lineage")}

{_disagreement_prose(disagree, missing)}

{caveat('''<p><b>What this cell typing does not establish.</b> Labels are transferred from a
reference built on different patients, so a population present in this tumour but absent from
the atlas cannot be named \u2014 it can only appear as low-confidence or unassigned. Cell
<em>states</em> (exhaustion, activation, EMT) are largely out of reach: the panel carries too few
of the genes that define them, and this analysis deliberately stops at lineage rather than
reporting state calls the data cannot support. Finally, the Xenium segmentation is nucleus
expansion, so transcripts from a tightly apposed neighbour are assigned to the wrong cell at
exactly the boundaries that matter most for the interface analysis that follows.</p>''')}
"""
    return body


def _disagreement_prose(disagree: pd.DataFrame, missing: pd.DataFrame) -> str:
    if len(disagree) == 0:
        return note("<p>Every cluster received the same broad call from both routes.</p>")
    rows = [[str(i), str(r["transferred_majority"]), str(r["marker_label"]),
             f"{int(r['n_cells']):,}", fmt(r["mean_transfer_confidence"], 2)]
            for i, r in disagree.iterrows()]
    miss = ""
    if len(missing):
        worst = missing.sort_values("n_found").iloc[0]
        miss = (f"<p>Panel coverage explains part of this: the <b>{worst['lineage']}</b> marker "
                f"set has only {int(worst['n_found'])} of {int(worst['n_requested'])} canonical "
                f"genes on the panel, so the marker route is working with less evidence than the "
                f"transfer route for that lineage.</p>")
    unass = disagree[disagree["transferred_majority"].astype(str) == "Unassigned"]
    extra = ""
    if len(unass):
        u = unass.iloc[0]
        extra = (f"<p>Cluster {u.name} is the most informative disagreement: the reference "
                 f"transfer leaves it <code>Unassigned</code> at a mean confidence of "
                 f"{fmt(u['mean_transfer_confidence'], 2)}, below the floor, while the marker "
                 f"route calls it <b>{u['marker_label']}</b>. That is the expected signature of "
                 f"a population the reference's coarse labels do not contain \u2014 the transfer "
                 f"route cannot name what it was never shown, and declining to guess is the "
                 f"correct behaviour. Here the second route recovers it. Neither route alone "
                 f"would have.</p>")
    return (f"{table(['Cluster', 'Reference transfer', 'Marker call', 'Cells', 'Mean confidence'], rows, numeric_from=3)}"
            f"{note(f'<p>These {len(disagree)} clusters carry a caveat, not a label. {miss}</p>{extra}', label='Flagged')}")


def page_spatial(d: dict) -> str:
    pairs, niche, iface = d["pairs"], d["niche"], d["iface"]
    seg = pairs.dropna(subset=["z"]).sort_values("z")
    most_seg, most_mix = seg.iloc[0], seg.iloc[-1]
    enr = iface[iface["cell_type"] != "Cancer Epithelial"].sort_values(
        "log2_enrichment", ascending=False)
    sig = enr[enr["q_bh"] < 0.05]
    irows = [[r["cell_type"], pct(r["frac_within_interface"], 1), pct(r["null_mean"], 1),
              fmt(r["log2_enrichment"]), fmt(r["z_score"], 1), _p(r["q_bh"])]
             for _, r in enr.iterrows()]
    nrows = [[str(i), ", ".join(f"{c} {v:.0%}" for c, v in
                                niche.loc[i].sort_values(ascending=False).head(3).items())]
             for i in niche.index]

    body = f"""
<p>Cell typing answers <em>what is here</em>. The tissue question is <em>what is next to what</em>,
and it needs a different set of tools: a spatial graph, a null model that respects the geometry,
and effect sizes rather than p-values alone. With {int(d['sc']['n_cells'].iloc[0]):,} cells,
almost any arrangement is "significant"; the number that matters is how far from random it is.</p>

<h2>The spatial graph</h2>

<p>Neighbours are defined by Delaunay triangulation with edges longer than
{int(CFG['spatial']['max_edge_um'])}&nbsp;\u00b5m removed. This choice matters more than
it looks. A fixed-<i>k</i> graph spans a few microns inside a dense tumour nest and tens of
microns in loose stroma, so a single <i>k</i> silently asks a different biological question in
each compartment. Delaunay adapts to local density; the distance cap stops cells on opposite
sides of a tissue void from being called neighbours.</p>

<h2>Which cell types touch each other?</h2>

{figure("assets/figures/fig07_neighborhood_enrichment.png",
        "<b>Neighbourhood enrichment.</b> Z-scores compare observed adjacency counts against "
        "a null that permutes cell-type labels while holding every cell position fixed, so the "
        "null preserves tissue density and shape exactly. Positive (blue) means the pair is "
        "found adjacent more often than chance; negative (red) means the pair is segregated.",
        alt="Heatmap and ranked bar chart of cell-type neighbourhood enrichment z-scores")}

<p>The strongest segregation is <b>{most_seg['pair']}</b> (z&nbsp;=&nbsp;{fmt(most_seg['z'], 0)})
and the strongest co-location is <b>{most_mix['pair']}</b> (z&nbsp;=&nbsp;{fmt(most_mix['z'], 0)}).
The diagonal of the matrix is strongly positive throughout, which is the expected signature of a
compartmentalised tissue: cells of a type are overwhelmingly adjacent to their own kind. In a
tumour that pattern is the quantitative statement of "nests and stroma", and its strength is the
thing worth comparing between patients.</p>

<h2>Recurrent niches</h2>

<p>Pairwise adjacency misses multi-way structure \u2014 a region where three populations are
interleaved is not described by any single pair. Niches capture it directly: for every cell,
take the cell-type composition of its {int(CFG['spatial']['niche_k'])} nearest neighbours,
and cluster those composition vectors across the whole section. Cells land in the same niche when
their local <em>mixture</em> is alike, regardless of where in the tissue they sit.</p>

{figure("assets/figures/fig08_niches.png",
        "<b>Niches and their composition.</b> (a) Each cell coloured by its niche; the niches "
        "reconstruct recognisable tissue architecture \u2014 tumour nest interiors, their "
        "boundaries, and distinct stromal and immune territories \u2014 without the algorithm "
        "ever being given a coordinate beyond the neighbour graph. (b) Mean neighbourhood "
        "composition of each niche, which is what makes the labels in (a) interpretable.",
        alt="Tissue map of spatial niches beside a stacked composition bar chart")}

{table(["Niche", "Dominant composition"], nrows, numeric_from=99)}

<h2>The tumour\u2013immune interface</h2>

<p>The clinically interesting question is not whether immune cells are present but whether they
reach the tumour. Each non-tumour cell is assigned its Euclidean distance to the nearest cell of
the tumour compartment, and the test asks whether a population is over- or under-represented
within {int(CFG['spatial']['interface_um'])}&nbsp;\u00b5m of that boundary relative to
a null that permutes cell-type labels across fixed positions
({int(CFG['spatial']['n_perms']):,} permutations, Benjamini\u2013Hochberg corrected).</p>

{figure("assets/figures/fig09_tumour_immune_interface.png",
        "<b>Distance to tumour and interface enrichment.</b> (a) Every cell coloured by distance "
        "to the nearest tumour cell, which makes the boundary geometry explicit. (b) Distance "
        "distributions per population. (c) Enrichment within the interface band against the "
        "permuted null, with BH-adjusted significance.",
        alt="Distance map, per-population distance distributions and interface enrichment bars")}

{table(["Population", "Observed in band", "Expected (null)", "log\u2082 enrichment", "z", "q (BH)"], irows)}

<p>{len(sig)} of {len(enr)} non-tumour populations differ significantly from the null at
q&nbsp;&lt;&nbsp;0.05. {_interface_prose(enr)}</p>

{caveat('''<p><b>Read these as description, not mechanism.</b> Proximity is not interaction: two
cells 20&nbsp;\u00b5m apart in a fixed section may never have signalled to each other, and
ligand\u2013receptor inference on a targeted panel would add precision without adding evidence.
The permutation null holds positions fixed and shuffles labels, so it controls for tissue
geometry and density but not for the possibility that cell typing errors are themselves
spatially structured \u2014 segmentation bleed at the tumour boundary is exactly such an error,
and it would act to inflate apparent immune presence at the interface. One section of one tumour
supports no claim about breast cancer in general.</p>''')}
"""
    return body


def _p(q: float) -> str:
    if q < 1e-3:
        return "&lt; 0.001"
    return f"{q:.3f}"


def _interface_prose(enr: pd.DataFrame) -> str:
    up = enr[(enr["log2_enrichment"] > 0) & (enr["q_bh"] < 0.05)]
    dn = enr[(enr["log2_enrichment"] < 0) & (enr["q_bh"] < 0.05)]
    bits = []
    if len(up):
        t = up.iloc[0]
        bits.append(f"<b>{t['cell_type']}</b> is the most strongly enriched population at the "
                    f"boundary (log<sub>2</sub>&nbsp;{fmt(t['log2_enrichment'])}, "
                    f"{pct(t['frac_within_interface'], 1)} of its cells within the band against "
                    f"{pct(t['null_mean'], 1)} expected)")
    if len(dn):
        t = dn.iloc[-1]
        bits.append(f"<b>{t['cell_type']}</b> is the most strongly excluded "
                    f"(log<sub>2</sub>&nbsp;{fmt(t['log2_enrichment'])})")
    if not bits:
        return "No population departs from the null after correction."
    return (" and ".join(bits) + ". The direction and magnitude together are the result; the "
            "p-value on this many cells is not informative on its own.")


def page_methods(d: dict, repo: str) -> str:
    thr = d["thr"]
    trows = [[k.replace("_", " "), f"{v:g}"] for k, v in thr.items()]
    body = f"""
<h2>Data</h2>
<p>Two public datasets, both downloaded by <code>scripts/download_data.sh</code>, which also
writes a <code>SHA256SUMS</code> manifest so a rerun can be verified byte-for-byte.</p>
<ul>
<li><b>10x Genomics Xenium In Situ, FFPE Human Breast Cancer, Replicates 1 and 2.</b> A
{int(d['sc']['n_panel_genes'].iloc[0])}-gene breast panel on two serial sections of one tumour
block; {int(d['sc']['n_cells'].sum()):,} segmented cells in total. Described in Janesick
et&nbsp;al., <i>Nature Communications</i> 2023.</li>
<li><b>Breast-cancer scRNA-seq atlas</b> of Wu et&nbsp;al., <i>Nature Genetics</i> 2021
(GEO accession GSE176078), used as the label-transfer reference.</li>
</ul>

<h2>Pipeline</h2>
<pre><code>bash scripts/download_data.sh     # public data + checksums
make qc                           # stage 1 -> results/tables/qc_scorecard.csv
make annotate                     # stage 2 -> data/processed/Rep1_annotated.h5ad
make spatial                      # stage 3 -> results/tables/immune_distance_stats.csv
python scripts/04_build_site.py   # this website
make test                         # unit tests, no data required</code></pre>

<h3>Stage 1 \u2014 quality control</h3>
<p>Counts, genes, control-feature fraction, cell and nucleus area and transcript density are
computed per cell from the raw matrix with negative-control features retained (they are dropped
only after the control rate has been computed from them). Run-level rates are recomputed from
the matrix rather than copied from the vendor <code>metrics_summary.csv</code>, which is carried
alongside for comparison. Regional artefacts are tested by binning the section into
{int(thr['grid_bin_um'])}&nbsp;\u00b5m squares and computing Moran's I with a queen-contiguity
weight matrix and a 199-permutation null.</p>

<h3>Stage 2 \u2014 cell typing</h3>
<p>Counts are normalised to a fixed total per cell and log1p-transformed; the full panel is used
as features, since variable-gene selection is meaningless on a targeted panel where every gene
was chosen for informativeness. PCA to {int(CFG['cluster']['n_pcs'])} components, Leiden
clustering over a resolution sweep with the chosen value recorded in the config, and Harmony
integration of Xenium with the reference before a distance-weighted k-nearest-neighbour vote.
Cells below the confidence floor are left unassigned. Marker scoring uses
<code>scanpy.tl.score_genes</code> against curated per-lineage sets, aggregated to clusters.</p>

<h3>Stage 3 \u2014 spatial analysis</h3>
<p>Delaunay graph with a distance cap; neighbourhood enrichment by label permutation with
positions held fixed; niches by k-means on fixed-<i>k</i> neighbourhood composition vectors;
interface enrichment by permutation with Benjamini\u2013Hochberg correction across populations.
Distances are Euclidean in microns on the vendor's cell centroids.</p>

<h2>Parameters</h2>
<p>The complete set, as applied. These are read from <code>config/config.yaml</code> at run time
and written to <code>results/tables/</code> with the outputs they produced.</p>
{table(["Parameter", "Value"], trows)}

<h2>Reproducibility</h2>
<ul>
<li><b>Pinned environment.</b> <code>environment.yml</code> for conda;
<code>pyproject.toml</code> for pip.</li>
<li><b>Deterministic.</b> One seed in the config threads through Leiden, Harmony, k-means and
every permutation test.</li>
<li><b>Tested.</b> {_test_count()} unit tests on synthetic fixtures with known answers, run in CI
on every push. The fixtures include a section with two spatially separated compartments, so
niche recovery and the interface null can be checked against a ground truth.</li>
<li><b>Checksummed inputs.</b> The download script writes SHA-256 sums for every file it
fetches.</li>
<li><b>No hand-edited numbers.</b> Every figure is produced by a script, and every number on
this website is injected from <code>results/tables/</code> when the site is built.</li>
</ul>

<h2>What this analysis cannot tell you</h2>
<ul>
<li><b>One tumour.</b> Two sections of one block. Nothing here generalises to breast cancer as a
disease; the value is in the method, not the biological claim.</li>
<li><b>No morphology channel.</b> H&amp;E-to-DAPI registration, and therefore any
pathologist-annotated region analysis, is outside what this public bundle supports.</li>
<li><b>Segmentation is nucleus expansion.</b> Transcript misassignment between tightly apposed
cells is unquantified here and is worst exactly at the tumour boundary.</li>
<li><b>Lineage, not state.</b> The panel cannot support exhaustion, activation or EMT calls, so
none are made.</li>
<li><b>Proximity is not interaction.</b> Every spatial result is a description of arrangement.</li>
<li><b>Assumption diagnostics were not run.</b> The permutation tests are distribution-free by
construction; no further residual or normality checks were performed.</li>
</ul>

<h2>References</h2>
<ol>
<li>Janesick A, et&nbsp;al. High resolution mapping of the tumor microenvironment using
integrated single-cell, spatial and in situ analysis. <i>Nature Communications</i> 14, 8353
(2023).</li>
<li>Wu SZ, et&nbsp;al. A single-cell and spatially resolved atlas of human breast cancers.
<i>Nature Genetics</i> 53, 1334\u20131347 (2021). GEO: GSE176078.</li>
<li>Palla G, et&nbsp;al. Squidpy: a scalable framework for spatial omics analysis.
<i>Nature Methods</i> 19, 171\u2013178 (2022).</li>
<li>Korsunsky I, et&nbsp;al. Fast, sensitive and accurate integration of single-cell data with
Harmony. <i>Nature Methods</i> 16, 1289\u20131296 (2019).</li>
<li>Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis.
<i>Genome Biology</i> 19, 15 (2018).</li>
</ol>

<p>Code: <a href="{repo}">{repo.replace('https://', '')}</a>. Author:
<a href="{AUTHOR_SITE}">{AUTHOR_SITE.replace('https://', '')}</a>.</p>
"""
    return body


def page_about(repo: str) -> str:
    return f"""
<p>This project is a self-contained demonstration of how I approach spatial transcriptomics
data: start from the public raw output rather than a processed object, make the quality criteria
explicit before looking at the result, run the annotation twice by independent routes, and state
plainly what the data cannot support.</p>

<h2>Why these three stages</h2>

<p>They are the stages where spatial analyses actually go wrong. QC on in situ data is
unforgiving \u2014 a background rate or a segmentation artefact propagates silently through
everything downstream, and by the time it surfaces it looks like biology. Cell typing on a
targeted panel is a judgement call that a pipeline cannot make for you: the panel sets a hard
ceiling on resolvable identity, and pretending otherwise produces confident labels on evidence
that does not exist. And spatial statistics is where a weak null model turns tissue geometry
into a false discovery; permuting labels while holding positions fixed is the minimum defensible
choice, and it is not the default in most published analyses.</p>

<h2>What is in the repository</h2>

<ul>
<li><code>src/xtme/</code> \u2014 a small library: readers, QC, annotation, spatial statistics
and a shared figure style. Every statistical function is importable and tested.</li>
<li><code>scripts/</code> \u2014 four numbered stages plus the download script. Each writes its
tables and figures and nothing else.</li>
<li><code>config/config.yaml</code> \u2014 every threshold and parameter in the analysis.</li>
<li><code>tests/</code> \u2014 unit tests on synthetic fixtures with known answers, run in CI on
every push without any data download.</li>
<li><code>results/</code> \u2014 the figures and tables behind every claim on this site,
committed so the repository is readable without running anything.</li>
</ul>

<h2>Reusing it</h2>

<p>The pipeline is parameterised on the config file, so pointing it at a different Xenium run
means editing paths and, if the tissue differs, swapping the marker sets in
<code>src/xtme/annotate.py</code> and the reference in the config. The QC and spatial modules are
tissue-agnostic.</p>

{note(f'''<p>Everything here is MIT-licensed and reproducible from a clean clone:
<a href="{repo}">{repo.replace("https://", "")}</a>. Issues and corrections are welcome \u2014
if something in the statistics looks wrong to you, I would rather know.</p>''')}

<h2>Contact</h2>

<p class="author-line">Ashok Dwivedi \u2014 <a href="{AUTHOR_SITE}">{AUTHOR_SITE.replace('https://', '')}</a></p>
<p>Background, other projects and contact details are on
<a href="{AUTHOR_SITE}">{AUTHOR_SITE.replace('https://', '')}</a>.</p>
"""


if __name__ == "__main__":
    main()
