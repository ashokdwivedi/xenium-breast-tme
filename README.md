# xenium-breast-tme

**Reproducible QC, reference-based cell typing and spatial neighbourhood analysis of a public
10x Genomics Xenium In Situ breast-cancer dataset.**

[![ci](https://github.com/ashokdwivedi/xenium-breast-tme/actions/workflows/ci.yml/badge.svg)](https://github.com/ashokdwivedi/xenium-breast-tme/actions/workflows/ci.yml)
[![license: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

**→ Full write-up with figures: https://ashokdwivedi.github.io/xenium-breast-tme/**
**→ Author: [ashokdwivedi.com](https://ashokdwivedi.com)**

---

Two replicate sections of one FFPE human breast tumour, taken from raw vendor output to
spatial statistics on the tumour–immune interface. The point of the repository is not the
biological result — it is one tumour — but the shape of the analysis: explicit quality
criteria, a cross-replicate control for batch effects, cell typing run twice by independent
routes with the disagreement reported, and spatial tests against a null model that respects
tissue geometry.

| | |
|---|---|
| **Data** | 10x Xenium FFPE Human Breast Cancer, Rep 1 + Rep 2 (313-gene panel, 286,532 cells) |
| **Reference** | Wu et al. 2021 breast-cancer scRNA-seq atlas (GEO GSE176078) |
| **Stack** | scanpy · squidpy · harmonypy · scikit-learn |
| **Tests** | 21 unit tests on synthetic fixtures, run in CI without any data download |

## Quickstart

```bash
git clone https://github.com/ashokdwivedi/xenium-breast-tme.git
cd xenium-breast-tme
conda env create -f environment.yml && conda activate xenium-breast-tme
pip install -e .

make test          # unit tests — no data needed, ~90 s
bash scripts/download_data.sh    # ~1.4 GB of public data + SHA256 manifest
make all           # QC → cell typing → spatial analysis
python scripts/04_build_site.py  # rebuild the website in docs/
```

## What each stage does

### Stage 1 — quality control (`make qc`)

Three levels, kept separate because they fail for different reasons:

- **Per cell.** Transcripts, genes, control-feature fraction, cell and nucleus area,
  transcript density. Every threshold lives in `config/config.yaml` and is written to
  `results/tables/qc_thresholds.csv` next to the verdict it produced.
- **Per run.** Negative-control probe and codeword rates recomputed from the count matrix
  rather than copied from the vendor summary (which is carried alongside for comparison),
  median depth, and cell yield — combined into a single `PASS`/`REVIEW` verdict.
- **Position-aware.** The section is binned into 100 µm squares and Moran's I is computed for
  each technical metric on a queen-contiguity grid. This is the check that separates a tissue
  gradient from a defocused tile or a decoding dropout, and it is the one most often skipped.

Both replicates pass. Negative-control probes account for 0.04–0.09% of panel signal.

**Cross-replicate concordance** is treated as a first-class result: pseudobulk expression over
the 313 panel genes correlates at r = 0.992 between sections, and the handful of genes that
differ more than two-fold are adipose and muscle markers — a sectioning difference, not a run
difference. That number is the floor on what a single-section fold-change can mean.

### Stage 2 — cell typing (`make annotate`)

Two independent routes, compared rather than merged:

- **Reference transfer.** Harmony-integrated k-NN vote from the Wu et al. atlas restricted to
  panel genes, with each reference label capped so an abundant population cannot win on prior
  alone. Cells below a confidence floor are labelled `Unassigned` rather than forced.
- **Marker scoring.** Leiden clusters scored against curated canonical marker sets; clusters
  whose best and second-best scores are close are flagged `ambiguous`.

Where the routes disagree is published, not smoothed over, along with per-lineage marker
coverage — a lineage the panel barely covers cannot be resolved by either route, and saying so
is part of the result.

### Stage 3 — spatial structure (`make spatial`)

- **Neighbourhood enrichment** on a Delaunay graph with a 60 µm edge cap. Fixed-*k* graphs span
  different physical distances in dense tumour versus loose stroma, which silently changes the
  question being asked; Delaunay adapts to local density.
- **Niches** by k-means on fixed-*k* neighbourhood composition vectors, so recurrent multi-way
  mixtures are recovered rather than only pairwise adjacency.
- **Tumour–immune interface** enrichment within 30 µm of the tumour compartment, tested against
  a label-permutation null that holds every cell position fixed, reported with effect sizes and
  BH-adjusted p-values.

### Stage 4 — website (`python scripts/04_build_site.py`)

Builds `docs/` for GitHub Pages. Every number in the prose is injected from
`results/tables/` at build time, so the write-up cannot drift from the analysis.

## Layout

```
src/xtme/          importable, tested library
  io.py            Xenium outs reader (keeps control features), reference streamer
  qc.py            per-cell / per-run / grid QC, Moran's I, pseudobulk
  annotate.py      normalisation, clustering, Harmony transfer, marker scoring
  spatial.py       graphs, neighbourhood enrichment, niches, interface statistics
  plotting.py      one figure style, colour bound once to cell-type identity
scripts/           01_qc · 02_annotate · 03_spatial · 04_build_site · download_data.sh
config/config.yaml every threshold and parameter in the analysis
tests/             unit tests on synthetic fixtures with known answers
results/           figures and tables, committed so the repo reads without running
docs/              the generated website
```

## Design decisions worth arguing about

- **Control features are kept at read time**, not dropped. Negative-control probe and codeword
  rates are the backbone of run-level QC; discarding them at load makes the most important
  quality metric unrecoverable.
- **Run metrics are recomputed, not trusted.** The vendor summary is carried alongside so the
  two can be compared, which is how you notice a bundle that has been subset or re-processed.
- **No variable-gene selection.** Every gene on a targeted panel was chosen for informativeness;
  HVG selection on 313 such genes discards signal to no benefit.
- **Harmony before the k-NN vote.** The dominant axis of variation between scRNA-seq and in situ
  data is the assay. Without integration, a nearest-neighbour vote recovers the batch.
- **Label permutation with fixed positions** as the spatial null. It preserves tissue density and
  shape exactly, so enrichment cannot be manufactured by geometry.
- **Effect sizes over p-values.** With 163,672 cells, nearly any arrangement is "significant".

## Limitations

One tumour, two sections — nothing here generalises to the disease. No morphology channel, so
H&E-to-DAPI registration and pathologist-annotated regions are out of scope for this public
bundle. Segmentation is nucleus expansion, so transcript misassignment is worst exactly at the
tumour boundary that the interface analysis depends on. The panel cannot support cell-state
calls (exhaustion, activation, EMT), so none are made. Proximity is not interaction.

## Citation

Data: Janesick et al., *Nat Commun* 14:8353 (2023); Wu et al., *Nat Genet* 53:1334 (2021),
GEO GSE176078. This repository is an independent methods demonstration and is not affiliated
with or endorsed by 10x Genomics.

## Licence

MIT — see [LICENSE](LICENSE).

---

Built by **Ashok Dwivedi** · [ashokdwivedi.com](https://ashokdwivedi.com)
