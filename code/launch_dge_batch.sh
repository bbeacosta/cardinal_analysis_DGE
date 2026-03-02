#!/usr/bin/env bash
set -euo pipefail

timestamp=$(date +%Y%m%d_%H%M%S)

############################
# 1) SETTINGS #
############################

# DNAnexus project & output parent folder
DX_PROJECT="project-Gx25k98J08pkk84J3V1JPPGY" 
OUT_PARENT="/Beatrice/results_CT3/${timestamp}"

# Where your repo is *on the machine you're launching from* (RAP RStudio session)
REPO_LOCAL="$HOME/analysis_code/cardinal_analysis_DGE"

# Entry scripts you will upload (from local machine to DNAnexus)
RUNSH_LOCAL="${REPO_LOCAL}/run.sh"
RWRAP_LOCAL="${REPO_LOCAL}/dge_run.R"

# Your main analysis script inside repo (will be *present in job* after we upload the repo tarball)
MAIN_R_SCRIPT="code/scDEA_LM_CT3_UKB_batch.R"  

# Main inputs stored on DNAnexus (use the DNAnexus paths, not local)
# (adjust these to where you keep them in the project)
FILTERED_COUNTS_DX="/Beatrice/results_CT3/filtered_counts_list.rds"
PCA_LIST_DX="/Beatrice/results_CT3/PCA_CT3/pca_list.rds"
MASTER_META_DX="/data/freeze3/F3_UKB_adata_obs_with_metadata.csv"
PHENOS_DX="/Daniela/data_F3/disease/cases_controls_all.tsv"

# Batch unit: cell types to run
# Option 1: explicit list
CELLTYPES=("B_CD5" "B_memory_IGHMlow" "HSC_MPP")
# Option 2 (recommended later): generate this list from filtered_counts_list.rds inside R, and pass --celltype from launcher.

# Job resources
INSTANCE_TYPE="mem3_ssd1_v2_x48"
PRIORITY="high"

##########################################
# 2) PREPARE OUTPUT + UPLOAD CODE BUNDLE #
##########################################

dx mkdir -p "${DX_PROJECT}:${OUT_PARENT}"

# Bundle your repo into a tar.gz so the job gets the exact code version you have
# (this avoids needing git clone inside the job)
echo "Bundling repo from ${REPO_LOCAL}"
tar -czf cardinal_analysis_DGE_${timestamp}.tar.gz -C "$(dirname "$REPO_LOCAL")" "$(basename "$REPO_LOCAL")"

# Upload tarball + run.sh + dge_run.R
CODEDIR="${OUT_PARENT}/code"
dx mkdir -p "${DX_PROJECT}:${CODEDIR}"

echo "Uploading code bundle + wrappers..."
dx upload "cardinal_analysis_DGE_${timestamp}.tar.gz" --path "${DX_PROJECT}:${CODEDIR}/cardinal_analysis_DGE.tar.gz" --brief
dx upload "${RUNSH_LOCAL}" --path "${DX_PROJECT}:${CODEDIR}/run.sh" --brief
dx upload "${RWRAP_LOCAL}" --path "${DX_PROJECT}:${CODEDIR}/dge_run.R" --brief

#####################################
# 3) SUBMIT ONE JOB PER CELL TYPE   #
#####################################

for ct in "${CELLTYPES[@]}"; do
  OUTDIR="${OUT_PARENT}/${ct}"
  dx mkdir -p "${DX_PROJECT}:${OUTDIR}"

  echo "Submitting celltype job: ${ct}"

  dx run swiss-army-knife \
    --brief \
    --destination "${DX_PROJECT}:${OUTDIR}" \
    --instance-type "${INSTANCE_TYPE}" \
    --priority "${PRIORITY}" \
    --name "DGE_${ct}_${timestamp}" \
    -iin "${DX_PROJECT}:${CODEDIR}/cardinal_analysis_DGE.tar.gz" \
    -iin "${DX_PROJECT}:${CODEDIR}/run.sh" \
    -iin "${DX_PROJECT}:${CODEDIR}/dge_run.R" \
    -iin "${DX_PROJECT}:${FILTERED_COUNTS_DX}" \
    -iin "${DX_PROJECT}:${PCA_LIST_DX}" \
    -iin "${DX_PROJECT}:${MASTER_META_DX}" \
    -iin "${DX_PROJECT}:${PHENOS_DX}" \
    -icmd "bash run.sh \
      --celltype '${ct}' \
      --repo_tar 'cardinal_analysis_DGE.tar.gz' \
      --main_r '${MAIN_R_SCRIPT}' \
      --filtered_counts '$(basename ${FILTERED_COUNTS_DX})' \
      --pca_list '$(basename ${PCA_LIST_DX})' \
      --master_meta '$(basename ${MASTER_META_DX})' \
      --phenos '$(basename ${PHENOS_DX})'" \
    --yes
done

echo "All jobs submitted. Outputs under: ${OUT_PARENT}"
