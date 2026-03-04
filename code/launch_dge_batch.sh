#!/usr/bin/env bash
set -euo pipefail

############################################
# USER SETTINGS (EDIT THESE ONLY IF NEEDED)
############################################

# DNAnexus destination where outputs will be uploaded
DX_OUT="Beatrice/results_CT3/results_DGE_batch_fullrun"

# Where to store the code bundle on DNAnexus
DX_CODEDIR="Beatrice/dge_code_bundle"

# DNAnexus inputs (must exist in your project)
DX_FILTERED_COUNTS="Beatrice/results_CT3/filtered_counts_list.rds"
DX_PCA_LIST="Beatrice/results_CT3/PCA_CT3/pca_list.rds"
DX_MASTER_META="data/freeze3/F3_UKB_adata_obs_with_metadata.csv"
DX_PHENOS="Daniela/data_F3/disease/cases_controls_all.tsv"

# Your main pipeline script path INSIDE the repo (after unpack)
MAIN_R="code/scDEA_LM_CT3_UKB_batch.R"   

# Instance type for the job
INSTANCE="mem3_ssd1_v2_x48"

# Job name prefix
STAMP="$(date +%Y%m%d_%H%M%S)"
JOB_NAME="DGE_CT3_fullrun_${STAMP}"

############################################
# BUILD REPO TAR (from your current repo)
############################################

# Expect to be run from inside the repo folder
# i.e. ~/analysis_code/cardinal_analysis_DGE
REPO_NAME="cardinal_analysis_DGE"
TAR_NAME="${REPO_NAME}.tar.gz"

echo "==> Creating repo tarball: ${TAR_NAME}"
cd ..
tar -czf "${TAR_NAME}" "${REPO_NAME}"
cd "${REPO_NAME}"

############################################
# UPLOAD BUNDLE + SCRIPTS
############################################

echo "==> Creating output + code dirs (mkdir is safe; won't delete existing)"
dx mkdir -p "${DX_CODEDIR}"
dx mkdir -p "${DX_OUT}"

echo "==> Uploading scripts to ${DX_CODEDIR}"
dx upload launch_dge_batch.sh --path "${DX_CODEDIR}/launch_dge_batch.sh" --brief
dx upload run.sh              --path "${DX_CODEDIR}/run.sh"              --brief
dx upload dge_run.R           --path "${DX_CODEDIR}/dge_run.R"           --brief

echo "==> Uploading repo tarball to ${DX_CODEDIR}"
dx upload "../${TAR_NAME}" --path "${DX_CODEDIR}/${TAR_NAME}" --brief

############################################
# LAUNCH ONE JOB (FULL RUN)
############################################

echo "==> Launching ONE job: ${JOB_NAME}"

dx run --brief --destination "${DX_OUT}" swiss-army-knife \
  -iin="${DX_CODEDIR}/run.sh" \
  -iin="${DX_CODEDIR}/dge_run.R" \
  -iin="${DX_CODEDIR}/${TAR_NAME}" \
  -iin="${DX_FILTERED_COUNTS}" \
  -iin="${DX_PCA_LIST}" \
  -iin="${DX_MASTER_META}" \
  -iin="${DX_PHENOS}" \
  --instance-type "${INSTANCE}" \
  --name "${JOB_NAME}" \
  -icmd="bash run.sh \
    --repo_tar ${TAR_NAME} \
    --main_r ${MAIN_R} \
    --filtered_counts $(basename ${DX_FILTERED_COUNTS}) \
    --pca_list $(basename ${DX_PCA_LIST}) \
    --master_meta $(basename ${DX_MASTER_META}) \
    --phenos $(basename ${DX_PHENOS}) \
    --dx_out ${DX_OUT}" \
  --priority high \
  --yes

echo "==> Submitted. Check Monitor tab in DNAnexus for job: ${JOB_NAME}"
