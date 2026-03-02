#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
run.sh --celltype CT --repo_tar cardinal_analysis_DGE.tar.gz --main_r code/script.R
       --filtered_counts filtered_counts_list.rds --pca_list pca_list.rds
       --master_meta F3_UKB_adata_obs_with_metadata.csv --phenos cases_controls_all.tsv
EOF
  exit 1
}

CELLTYPE=""
REPO_TAR=""
MAIN_R=""
FILTERED=""
PCA=""
META=""
PHENOS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --celltype) CELLTYPE="$2"; shift 2;;
    --repo_tar) REPO_TAR="$2"; shift 2;;
    --main_r) MAIN_R="$2"; shift 2;;
    --filtered_counts) FILTERED="$2"; shift 2;;
    --pca_list) PCA="$2"; shift 2;;
    --master_meta) META="$2"; shift 2;;
    --phenos) PHENOS="$2"; shift 2;;
    *) echo "Unknown argument: $1"; usage;;
  esac
done

[[ -n "$CELLTYPE" && -n "$REPO_TAR" && -n "$MAIN_R" ]] || usage

echo "Working directory: $(pwd)"
echo "Files present:"
ls -lah

# unpack repo
tar -xzf "$REPO_TAR"

# repo folder name after unpack (must match tar creation)
cd cardinal_analysis_DGE

echo "Repo unpacked. Current dir: $(pwd)"
ls -lah

# run R wrapper
Rscript dge_run.R \
  --celltype "${CELLTYPE}" \
  --main_r "${MAIN_R}" \
  --filtered_counts "../${FILTERED}" \
  --pca_list "../${PCA}" \
  --master_meta "../${META}" \
  --phenos "../${PHENOS}"

echo "Job finished for celltype: ${CELLTYPE}"
