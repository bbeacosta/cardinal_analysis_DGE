#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
run.sh --repo_tar cardinal_analysis_DGE.tar.gz --main_r code/script.R
       --filtered_counts filtered_counts_list.rds --pca_list pca_list.rds
       --master_meta F3_UKB_adata_obs_with_metadata.csv --phenos cases_controls_all.tsv
       [--dx_out Beatrice/results_CT3/results_DGE_batch_fullrun]
EOF
  exit 1
}

REPO_TAR=""
MAIN_R=""
FILTERED=""
PCA=""
META=""
PHENOS=""
DX_OUT=""   # optional DNAnexus destination folder (project-relative)

while [[ $# -gt 0 ]]; do
  case "$1" in
    --repo_tar) REPO_TAR="$2"; shift 2;;
    --main_r) MAIN_R="$2"; shift 2;;
    --filtered_counts) FILTERED="$2"; shift 2;;
    --pca_list) PCA="$2"; shift 2;;
    --master_meta) META="$2"; shift 2;;
    --phenos) PHENOS="$2"; shift 2;;
    --dx_out) DX_OUT="$2"; shift 2;;
    *) echo "Unknown argument: $1"; usage;;
  esac
done

[[ -n "$REPO_TAR" && -n "$MAIN_R" && -n "$FILTERED" && -n "$PCA" && -n "$META" && -n "$PHENOS" ]] || usage

echo "Working directory: $(pwd)"
echo "Files present:"
ls -lah

# DNAnexus upload directory (only this gets uploaded automatically)
OUT_ROOT="/home/dnanexus/out"
mkdir -p "${OUT_ROOT}"

# unpack repo
tar -xzf "$REPO_TAR"
cd cardinal_analysis_DGE

echo "Repo unpacked. Current dir: $(pwd)"
ls -lah

# IMPORTANT: force your pipeline to write into /home/dnanexus/out
# Your R script currently uses base_dir <- "/home/rstudio-server"
# Override with an env var your R script can read (recommended), OR
# simplest: create a symlink so "/home/rstudio-server" points into out.
mkdir -p /home/rstudio-server
ln -sfn "${OUT_ROOT}" /home/rstudio-server

# run R wrapper
Rscript dge_run.R \
  --main_r "${MAIN_R}" \
  --filtered_counts "../${FILTERED}" \
  --pca_list "../${PCA}" \
  --master_meta "../${META}" \
  --phenos "../${PHENOS}"

echo "R finished. Contents of ${OUT_ROOT}:"
find "${OUT_ROOT}" -maxdepth 3 -type d -print
find "${OUT_ROOT}" -maxdepth 3 -type f -print

# If user provided a DNAnexus destination folder, upload everything there.
# (If not, DNAnexus will still upload to the job --destination automatically.)
if [[ -n "${DX_OUT}" ]]; then
  echo "Uploading all outputs to: ${DX_OUT}"
  dx-upload-all-outputs --destination "${DX_OUT}" --recursive
fi

echo "Job finished"
