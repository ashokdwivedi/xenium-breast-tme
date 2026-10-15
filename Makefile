# ---------------------------------------------------------------------------
# xenium-breast-tme - end to end pipeline
#
#   make data      download the public Xenium runs and the scRNA-seq reference
#   make qc        stage 1: per-cell, per-run and spatial quality control
#   make annotate  stage 2: clustering + reference label transfer + markers
#   make spatial   stage 3: neighbourhoods, niches, tumour-immune interface
#   make all       everything, in order
#   make test      unit tests on synthetic fixtures (no data required)
# ---------------------------------------------------------------------------
PYTHON ?= python
RAW     := data/raw
INTERIM := data/interim
PROC    := data/processed
FIG     := results/figures
TAB     := results/tables

REF_DIR     := $(RAW)/GSE176078/Wu_etal_2021_BRCA_scRNASeq
REP1_MATRIX := $(RAW)/Xenium_FFPE_Human_Breast_Cancer_Rep1/cell_feature_matrix/matrix.mtx.gz
QC_STAMP    := $(TAB)/qc_scorecard.csv
ANN_STAMP   := $(PROC)/Rep1_annotated.h5ad
SPA_STAMP   := $(PROC)/Rep1_spatial.h5ad

.PHONY: all data qc annotate spatial test lint clean clean-results docs-serve
.DEFAULT_GOAL := all

all: spatial

data: $(REP1_MATRIX)
$(REP1_MATRIX):
	bash scripts/download_data.sh $(RAW)

qc: $(QC_STAMP)
$(QC_STAMP): $(REP1_MATRIX) scripts/01_qc.py src/xtme/qc.py src/xtme/io.py config/config.yaml
	$(PYTHON) scripts/01_qc.py

annotate: $(ANN_STAMP)
$(ANN_STAMP): $(QC_STAMP) scripts/02_annotate.py src/xtme/annotate.py config/config.yaml
	$(PYTHON) scripts/02_annotate.py

spatial: $(SPA_STAMP)
$(SPA_STAMP): $(ANN_STAMP) scripts/03_spatial.py src/xtme/spatial.py config/config.yaml
	$(PYTHON) scripts/03_spatial.py

test:
	$(PYTHON) -m pytest -q

lint:
	$(PYTHON) -m ruff check src scripts tests

clean-results:
	rm -rf $(FIG)/* $(TAB)/*

clean: clean-results
	rm -rf $(INTERIM)/* $(PROC)/*
