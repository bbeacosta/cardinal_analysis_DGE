#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
run.sh --repo_tar cardinal_analysis_DGE.tar.gz --main_r code/script.R
       --filtered_counts filtered_counts_list.rds --pca_list pca_list.rds
       --master_meta F3_UKB_adata_obs_with_metadata.csv --phenos cases_controls_all.tsv
EOF
  exit 1
}

REPO_TAR=""
MAIN_R=""
FILTERED=""
PCA=""
META=""
PHENOS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --repo_tar) REPO_TAR="$2"; shift 2;;
    --main_r) MAIN_R="$2"; shift 2;;
    --filtered_counts) FILTERED="$2"; shift 2;;
    --pca_list) PCA="$2"; shift 2;;
    --master_meta) META="$2"; shift 2;;
    --phenos) PHENOS="$2"; shift 2;;
    *) echo "Unknown argument: $1"; usage;;
  esac
done

[[ -n "$REPO_TAR" && -n "$MAIN_R" && -n "$FILTERED" && -n "$PCA" && -n "$META" && -n "$PHENOS" ]] || usage

echo "Working directory: $(pwd)"
echo "Files present in cwd:"
ls -lah

# Ensure your script's hardcoded path exists
mkdir -p /home/rstudio-server

echo "Unpacking repo tar: $REPO_TAR"
tar -xzf "$REPO_TAR"

cd cardinal_analysis_DGE
echo "Repo unpacked. Current dir: $(pwd)"
ls -lah

echo "Running wrapper (wrapper is staged one dir up)"
Rscript ../dge_run.R \
  --celltype "ALL" \
  --main_r "${MAIN_R}" \
  --filtered_counts "../${FILTERED}" \
  --pca_list "../${PCA}" \
  --master_meta "../${META}" \
  --phenos "../${PHENOS}"

echo "R finished. Listing /home/rstudio-server:"
ls -lah /home/rstudio-server || true

#############################################
# DNAnexus upload: copy outputs to /home/dnanexus/out
# Anything under /home/dnanexus/out is auto-uploaded
#############################################
mkdir -p /home/dnanexus/out

# Copy only if they exist (avoid failing if a folder wasn't created)
for d in results_DGE results_DGE_csv results_DGE_csv_annotated; do
  if [ -d "/home/rstudio-server/${d}" ]; then
    echo "Copying /home/rstudio-server/${d} -> /home/dnanexus/out/${d}"
    cp -a "/home/rstudio-server/${d}" "/home/dnanexus/out/${d}"
  else
    echo "Not found (skipping): /home/rstudio-server/${d}"
  fi
done

echo "Final contents of /home/dnanexus/out:"
find /home/dnanexus/out -maxdepth 3 -type f | head -n 200

echo "Job finished."
