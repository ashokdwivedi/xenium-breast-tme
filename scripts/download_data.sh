#!/usr/bin/env bash
# Download all public inputs for the xenium-breast-tme analysis.
#
#   bash scripts/download_data.sh [DATA_DIR]
#
# Inputs (all public, no credentials required):
#   1. 10x Genomics Xenium In Situ, FFPE Human Breast Cancer, Replicates 1 and 2
#      Xenium Analyzer v1.0.1, 313-plex Human Breast panel.
#      https://www.10xgenomics.com/datasets  (CC BY 4.0)
#   2. Wu et al. 2021, Nat Genet - breast cancer scRNA-seq atlas, GEO GSE176078.
#      Used as the cell-type reference for label transfer.
set -euo pipefail

DATA_DIR="${1:-data/raw}"
XEN_BASE="https://cf.10xgenomics.com/samples/xenium/1.0.1"
GEO_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz"

mkdir -p "$DATA_DIR"

fetch () {  # fetch <url> <dest>
    local url="$1" dest="$2"
    if [[ -s "$dest" ]]; then
        echo "  [skip] $(basename "$dest") already present"
        return
    fi
    echo "  [get ] $(basename "$dest")"
    curl -fSL --retry 3 --retry-delay 5 -o "${dest}.part" "$url"
    mv "${dest}.part" "$dest"
}

echo "== Xenium FFPE Human Breast Cancer (Rep1, Rep2)"
for REP in Rep1 Rep2; do
    SAMPLE="Xenium_FFPE_Human_Breast_Cancer_${REP}"
    OUT="${DATA_DIR}/${SAMPLE}"
    mkdir -p "$OUT"
    for F in cell_feature_matrix.tar.gz cells.csv.gz metrics_summary.csv gene_panel.json analysis.tar.gz; do
        fetch "${XEN_BASE}/${SAMPLE}/${SAMPLE}_${F}" "${OUT}/${F}"
    done
    # unpack the matrix (barcodes/features/matrix.mtx) and the 10x secondary analysis
    [[ -d "${OUT}/cell_feature_matrix" ]] || tar -xzf "${OUT}/cell_feature_matrix.tar.gz" -C "$OUT"
    [[ -d "${OUT}/analysis" ]]            || tar -xzf "${OUT}/analysis.tar.gz"            -C "$OUT"
done

echo "== Reference scRNA-seq atlas (Wu et al. 2021, GSE176078)"
REF="${DATA_DIR}/GSE176078"
mkdir -p "$REF"
fetch "$GEO_URL" "${REF}/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz"
if [[ ! -d "${REF}/Wu_etal_2021_BRCA_scRNASeq" ]]; then
    tar -xzf "${REF}/GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz" -C "$REF"
fi

echo "== Checksums -> ${DATA_DIR}/SHA256SUMS"
( cd "$DATA_DIR" && find . -type f \( -name '*.gz' -o -name '*.json' -o -name '*.csv' \) \
    -not -name 'SHA256SUMS' | sort | xargs shasum -a 256 ) > "${DATA_DIR}/SHA256SUMS"

echo "Done. Tree:"
du -sh "${DATA_DIR}"/*
