#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
run.sh --repo_tar cardinal_analysis_DGE.tar.gz --main_r code/script.R
       --filtered_counts filtered_counts_list.rds --pca_list pca_list.rds
       --master_meta F3_UKB_adata_obs_with_metadata.csv --phenos cases_controls_all.tsv
       --dx_out Beatrice/results_CT3/results_DGE_batch_fullrun
EOF
  exit 1
}

REPO_TAR=""
MAIN_R=""
FILTERED=""
PCA=""
META=""
PHENOS=""
DX_OUT=""

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

[[ -n "$REPO_TAR" && -n "$MAIN_R" && -n "$FILTERED" && -n "$PCA" && -n "$META" && -n "$PHENOS" && -n "$DX_OUT" ]] || usage

echo "==> Working directory: $(pwd)"
echo "==> Inputs present:"
ls -lah

############################################
# Stage files exactly where your R script expects them
############################################

echo "==> Staging inputs under /home/rstudio-server (to match your script)"

mkdir -p /home/rstudio-server
cp -v "$FILTERED" /home/rstudio-server/filtered_counts_list.rds

mkdir -p /home/rstudio-server/PCA_CT3
cp -v "$PCA" /home/rstudio-server/PCA_CT3/pca_list.rds

cp -v "$META" /home/rstudio-server/F3_UKB_adata_obs_with_metadata.csv
cp -v "$PHENOS" /home/rstudio-server/disease_cases_controls_all.tsv

echo "==> Confirm staged files:"
ls -lah /home/rstudio-server
ls -lah /home/rstudio-server/PCA_CT3

############################################
# Unpack repo
############################################

echo "==> Unpacking repo tarball: ${REPO_TAR}"
tar -xzf "$REPO_TAR"

cd cardinal_analysis_DGE
echo "==> Repo unpacked: $(pwd)"
ls -lah

############################################
# Run R (wrapper)
############################################

echo "==> Running R wrapper (full run over all CTs in your script)"
Rscript dge_run.R --main_r "${MAIN_R}"

############################################
# Upload outputs to DNAnexus
############################################

echo "==> Uploading outputs from /home/rstudio-server to DNAnexus: ${DX_OUT}"

# Upload the main result folders if they exist
for d in /home/rstudio-server/results_DGE /home/rstudio-server/results_DGE_csv /home/rstudio-server/results_DGE_csv_annotated; do
  if [[ -d "$d" ]]; then
    echo "   -> Uploading folder: $d"
    dx upload -r "$d" --path "${DX_OUT}/$(basename "$d")" --brief
  else
    echo "   -> Not found (skipping): $d"
  fi
done

# Upload any extra PDFs you might have created elsewhere (optional)
# Example: PCA_CT3 outputs or heatmaps saved outside results_DGE*
# Uncomment if needed:
if [[ -d /home/rstudio-server/PCA_CT3 ]]; then
dx upload -r /home/rstudio-server/PCA_CT3 --path "${DX_OUT}/PCA_CT3" --brief
fi

echo "==> Done. Outputs should now be in: ${DX_OUT}/"
