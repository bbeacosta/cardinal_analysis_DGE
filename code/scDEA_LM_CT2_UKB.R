
#----------------------------------------- FREEZE 3 -----------------------------------------------

##################################### LINEAR MODEL FOR DGE ###################################
# Load necessary libraries
.libPaths()

# Install packages
# Ensure BiocManager is available
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Full package list
pkgs <- c(
  "edgeR",
  "limma",
  "tidyr",
  "dplyr",
  "ggplot2",
  "EnhancedVolcano",
  "optparse",
  "stringr",
  "tibble",
  "pheatmap",
  "EnsDb.Hsapiens.v86",
  "AnnotationDbi"
)

# Install everything (BiocManager handles both CRAN + Bioconductor)
BiocManager::install(pkgs, ask = FALSE, update = FALSE)

# Load libraries
library(edgeR)
library(limma)
library(tidyr)
library(dplyr)
library(ggplot2)
library(EnhancedVolcano)
library(optparse)
library(stringr)
library(tibble)
library(pheatmap)
# library(lme4) # for mixed model only
library(EnsDb.Hsapiens.v86)
library(AnnotationDbi)



# Define base directory
base_dir <- "/home/rstudio-server"
dir.create(base_dir, showWarnings = FALSE)



################ Load all necessary files for analyses ##################
# load metadatafiles
# traits_macrocat <- read.csv("/genesandhealth/red/DanielaZanotti/data/disease/name_disease.csv", header = TRUE) # extended disease name for plotting and disease apparatus macrocategories

# phenotypes case-controls
phenotypes <-  read.table("/home/rstudio-server/cases_controls_all.tsv", header = TRUE, sep = "\t", check.names = FALSE)      
str(phenotypes)
head(phenotypes)

# # load PCA object
# pca_list <- readRDS("/home/ivm/LMM_DGE_results/pca_list.rds")
# 
# # load processed counts matrix
# filtered_counts_list <- readRDS("/home/ivm/Desktop/Beatrice_Costa/results_DEA/filtered_counts_list_fixed.rds")

# master table - donor level (each row is a cell)
master_table_F3 <- read.csv("/home/rstudio-server/F3_UKB_adata_obs_with_metadata.csv")      # Master table contains: age, sex, ethnicity, date at PBMC collection
str(master_table_F3)
head(master_table_F3)

# # obs - cell level (each row is a cell, not aggregated by donor like pseudobulks) - for ID conversion with phenotypes
# obs_clean <- read.table("/genesandhealth/red/DanielaZanotti/data/Freeze3/obs_gh-ContaminationFree-noDoublets-annotated-QCed-counts_CLEAN.tsv", header = TRUE, sep = "\t", check.names = FALSE)   # Vacutainer ID
# str(obs_clean)
# head(obs_clean)

############### Start processing the pseudobulk files #########################
##### load Freeze 3 pseudobulk data CT2 #####
# Create counts directories
dir_bangpaki <- "/home/rstudio-server/celltype_2/ukb-qced-cells/"

# extract everything before the first dot
all_files <- list.files(dir_bangpaki, full.names = FALSE)
celltypes <- unique(sub("\\..*", "", all_files))
exclude <- c(
  "gene_exp_variance",
  "make_parquet",
  "old_h5ad",
  "old_tables",
  "parquet_dataset_counts_long",
  "parquet_dataset_long",
  "update_sample_id",
  "update_sample_ids"
)

celltypes <- setdiff(celltypes, exclude)
print(celltypes) # 33 celltypes

# remove those that are not cell types

# Set up ncell filter: at least 10 cells per donor per celltype (based on tables in: *ncells_x_donor.tsv)
# List all counts files
counts_files <- list.files(dir_bangpaki, pattern = "\\.count\\.agg_sum\\.tsv$", full.names = TRUE)

# Output directory for filtered counts
out_dir <- file.path(dir_bangpaki, "filtered_counts")
dir.create(out_dir, showWarnings = FALSE)

# Run loop to filter counts
filtered_counts_list <- list()

for (cf in counts_files) {
  
  # 1) Extract celltype name from filename
  celltype <- str_replace(basename(cf), "\\.count\\.agg_sum\\.tsv$", "")
  message("\n==============================")
  message("Processing: ", celltype)
  
  # 2) Match ncells file
  ncells_file <- file.path(dir_bangpaki, paste0(celltype, ".ncells_x_donor.tsv"))
  if (!file.exists(ncells_file)) {
    warning("No ncells file found for ", celltype, " — skipping")
    next
  }
  
  # 3) Read counts (genes x donors)
  counts <- read.delim(cf, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, check.names = FALSE)
  gene_col   <- counts[[1]]
  counts_mat <- counts[,-1, drop = FALSE]
  rownames(counts_mat) <- gene_col
  
  # Normalize counts donor IDs: -- and __ → -
  colnames(counts_mat) <- colnames(counts_mat) |>
    str_replace_all("--", "-") |>
    str_replace_all("__", "-")
  
  # 4) Read ncells
  ncells <- read.delim(ncells_file, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, check.names = FALSE)
  
  # Pick a donor ID column robustly
  id_candidates <- intersect(
    c("unique.id", "unique_id", "donor_uid", "donor_id", "donor", "sample_id"),
    colnames(ncells)
  )
  
  if (length(id_candidates) == 0) {
    warning("No suitable donor ID column found in ncells for ", celltype,
            " — columns available: ", paste(colnames(ncells), collapse = ", "))
    next
  }
  
  id_col <- id_candidates[1]
  message("Using ncells ID column: ", id_col)
  
  # Normalize ncells IDs and filter N_cells
  ncells_keep <- ncells |>
    mutate(
      .donor_id_norm = .data[[id_col]] |>
        str_replace_all("--", "-") |>
        str_replace_all("__", "-")
    ) |>
    dplyr::filter(N_cells >= 10)
  
  message("ncells rows total: ", nrow(ncells),
          " ; ncells_keep (N_cells>=10): ", nrow(ncells_keep))
  
  # If all IDs are NA/empty after normalisation, skip
  if (nrow(ncells_keep) == 0 || all(is.na(ncells_keep$.donor_id_norm)) ||
      all(ncells_keep$.donor_id_norm == "")) {
    warning("For ", celltype,
            ": no usable donor IDs after filtering/normalization — skipping")
    next
  }
  
  # 5) Match donors between counts and ncells
  donors_keep <- intersect(colnames(counts_mat), ncells_keep$.donor_id_norm)
  
  if (length(donors_keep) == 0) {
    warning("No donors matched for ", celltype)
    message("Example counts donors: ",
            paste(head(colnames(counts_mat), 5), collapse = ", "))
    message("Example ncells donors (normalized): ",
            paste(head(na.omit(ncells_keep$.donor_id_norm), 5), collapse = ", "))
    next
  }
  
  counts_filt <- counts_mat[, donors_keep, drop = FALSE]
  message("Donors/columns kept: ", ncol(counts_filt), "/", ncol(counts_mat))
  
  # 6) Save filtered counts
  out_file <- file.path(out_dir, paste0(celltype, "_counts_filtered.tsv"))
  write.table(
    cbind(gene_id = rownames(counts_filt), counts_filt),
    file = out_file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  
  # 7) Store in list for R session use
  filtered_counts_list[[celltype]] <- counts_filt
}

# After loop: list of filtered counts per cell type
names(filtered_counts_list)
sapply(filtered_counts_list, ncol)

# Now drop celltypes with no donors >= 10 cells
filtered_counts_list <- Filter(function(x) ncol(x) > 0, filtered_counts_list)
sapply(filtered_counts_list, ncol)

# Save file for safe upload later
saveRDS(filtered_counts_list, file = "/home/rstudio-server/filtered_counts_list.rds")

############### Start processing the metadata #########################
### Match phenotypes IDs with counts IDs using obs_clean
# Create column in obs_clean that combines chromium_run_id and then _ donor_id columns to match pseudobulks
# Create a new column combining chromium_run_id and donor_id
# obs_clean <- obs_clean %>%
#   mutate(unique_id = paste(chromium_run_id, donor_id, sep = "_"))
# colnames(obs_clean)
# head(obs_clean$unique_id)
# 
# # Reduce obs_clean (cell-level metadata) to one row per donor (only keep unique ids for unique_id, state and vacutainer_id)
# obs_id_map <- obs_clean %>%
#   dplyr::select(vacutainer_id, unique_id, state) %>%
#   distinct(vacutainer_id, .keep_all = TRUE)  # one per donor
# 
# # Join with phenotypes (donor-level metadata) --> vacutainer_id (obs_clean) with CARDINAL_ID_sample1 (mastertable) columns
# phenotypes <- phenotypes %>%
#   left_join(obs_id_map, by = c("CARDINAL_ID_sample1" = "vacutainer_id"))

# Filtered counts list
readRDS(filtered_counts_list, file = "/home/rstudio-server/filtered_counts_list.rds")

# Remove NAs from $unique_id column
phenotypes <- phenotypes[!is.na(phenotypes$eid), ]

# Check I have only non-caucasian donors and no duplicates in metadata
# Metadata are cell-level, so I need to transform them to donor-level data
unique(master_table_F3$ancestry)

master_table_F3_donorL <- master_table_F3 %>%
  dplyr::select(
    eid,
    donor_uid_tpd,
    tranche_id,
    pool_id,
    sex,
    age,
    bmi,
    smoking_status_combined, 
    smoking_status_numeric
  ) %>%
  distinct()

any(duplicated(master_table_F3_donorL$eid)) # should return FALSE

master_table_F3_donorL <- master_table_F3_donorL %>%
  dplyr::mutate(
    donor_uid_tpd_norm = donor_uid_tpd %>%
      str_replace_all("\r", "") %>%   # just in case there are hidden CR chars
      str_trim() %>%
      str_replace_all("__", "-") %>%  # handle any __ cases too
      str_replace_all("--", "-")      # convert -- to single dash
  )



# ================= HELPER FUNCTIONS =================
cpm_edgeR <- function(x, log = FALSE, prior.count = 0.25) {
  edgeR::cpm(x, log = log, prior.count = prior.count)
}

# # Volcano plot function
plot_volcano <- function(
    tt,
    test_name,
    ct,
    gene_label_col = "Gene",
    topN = 10
){
  
  df <- tt
  
  # Check required columns
  required_cols <- c("logFC", "P.Value", "adj.P.Val")
  if (!all(required_cols %in% colnames(df))) {
    stop("topTable object missing one of: logFC, P.Value, adj.P.Val")
  }
  
  # If gene column missing, create one
  if (!(gene_label_col %in% colnames(df))) {
    df[[gene_label_col]] <- rownames(df)
  }
  
  # Volcano categories
  df$category <- dplyr::case_when(
    df$adj.P.Val < 0.05 & abs(df$logFC) >= 1 ~ "FDR < 0.05 & |LFC|",
    df$adj.P.Val < 0.05                     ~ "FDR < 0.05",
    abs(df$logFC) >= 1                      ~ "|LFC|",
    TRUE                                    ~ "Not significant"
  )
  
  colors <- c(
    "FDR < 0.05 & |LFC|" = "red",
    "FDR < 0.05 "        = "blue",
    "|LFC|"         = "orange",
    "Not significant"        = "grey70"
  )
  
  # Pick top 10 genes (by FDR)
  top_labels <- df[order(df$adj.P.Val), ][1:min(topN, nrow(df)), ]
  
  p <- ggplot(df, aes(x = logFC,
                      y = -log10(P.Value),
                      color = category)) +
    geom_point(alpha = 0.8, size = 1.2) +
    scale_color_manual(values = colors) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_text(
      data = top_labels,
      aes(label = !!sym(gene_label_col)),
      size = 3,
      vjust = -0.4,
      color = "black",
      check_overlap = TRUE
    ) +
    theme_classic(base_size = 12) +
    labs(
      title = paste("DGE by", test_name, "-", ct),
      x = "log2 Fold Change",
      y = "-log10(p-value)",
      color = "Significance"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
  
  return(p)
}



# Create list to later compute variance explained by each covariate, for case/control counts, qc metrics
variance_explained_list <- list()
case_counts <- list()
qc_metrics <- list()

# ================= MAIN LOOP =================
# # Remove problematic cell type
# counts_list[["HSC_MPP"]] <- NULL
# counts_list[["ILC"]] <- NULL
# counts_list[["NK_prolif"]] <- NULL 
# counts_list[["Plasmablast"]] <- NULL
# counts_list[["T_CD8_prolif"]] <- NULL
# 

# ct <- "B_exhausted"   

for (ct in names(filtered_counts_list)) {
  message("\n===== CELL TYPE: ", ct, " =====")
  
  # Iteratively create results folders for specific outputs
  # Create folder for this cell type
  ct_dir <- file.path(base_dir, ct)
  data_dir  <- file.path(ct_dir, "data")
  plots_dir <- file.path(ct_dir, "plots")
  
  dir.create(ct_dir,  showWarnings = FALSE)
  dir.create(data_dir,  showWarnings = FALSE)
  dir.create(plots_dir, showWarnings = FALSE)
  
  
  ### Number of PCs to include as technical covariates
  # number of genetic/technical PCs to use instead of pool/tranche
  pc_scores_df <- as.data.frame(pca_list[[ct]]$x) # extract PCA score for each cell type
  pc_scores_df$unique_id <- rownames(pc_scores_df) #create unique_id column from rownames
  rownames(pc_scores_df) <- NULL
  n_PC_use <- 4
  pc_scores_df_name <- "pc_scores_df"   # data.frame with columns: unique_id, PC1..PCn
  
  # Input data
  counts_list <- filtered_counts_list    # genes x samples
  
  bio_covs <- c("ancestry", "scaled_age", "sex", "scaled_BMI")
  pc_prefix <- "PC"      # columns in pc_scores_df
  
  out_dir <- "DGE_results_LM"
  dir.create(out_dir, showWarnings = FALSE)
  
  counts <- counts_list[[ct]]
  if (is.null(counts) || !is.matrix(counts) || ncol(counts) < 4) {
    message("Skipping ", ct, ": invalid counts or too few samples.")
    next
  }
  # Remove NAs from $pool_id column in mastertable
  meta_all <- master_table_F3 [!is.na(master_table_F3$pool_id), ]
  
  # Combine pool_id and donor_id into a single column with an underscore
  meta_all$unique_id <- paste(meta_all$pool_id, meta_all$donor_id, sep = "_")
  
  # Check the first few to confirm
  head(meta_all$unique_id)
  
  # Align metadata
  shared_ids <- intersect(colnames(counts), meta_all$unique_id)
  counts <- counts[, shared_ids, drop = FALSE]
  meta_sub <- meta_all[match(shared_ids, meta_all$unique_id), , drop = FALSE]
  stopifnot(all(colnames(counts) == meta_sub$unique_id))
  
  # Make ancestry and sex as factors so that they are treated by limma as binary variables and not continuous variables
  meta_sub$sex <- factor(meta_sub$sex, levels = c(1, 2), labels = c("F", "M"))
  meta_sub$ancestry<- factor(meta_sub$ancestry)
  
  # Optional but recommended: set biologically meaningful reference
  # adjust levels as appropriate for your encoding
  levels(meta_sub$sex)
  meta_sub$sex <- relevel(meta_sub$sex, ref = "M")
  
  
  # Scale continuous variables so that they all have mean = 0 and sd = 1, to improve model stability and make effect sizes comparable
  # When covariates have different numeric scales, scaling makes all covariates comparable in scale, helping model algorithm converge faster and more reliably
  meta_sub$scaled_age <- as.numeric(scale(meta_sub$age_at_recruitment_first))
  meta_sub$scaled_BMI <- as.numeric(scale(meta_sub$BMI))
  
  
  # ----- Attach PCs as technical covariates -----
  # Attach PCs
  if (exists(pc_scores_df_name, envir = .GlobalEnv)) {
    pc_scores_df <- get(pc_scores_df_name, envir = .GlobalEnv)
    pc_cols <- grep(paste0("^", pc_prefix), colnames(pc_scores_df), value = TRUE)
    if (length(pc_cols) > 0) {
      pc_cols_use <- head(pc_cols, n_PC_use)
      pc_scores_sub <- pc_scores_df[match(meta_sub$unique_id, pc_scores_df$unique_id), pc_cols_use, drop = FALSE]
      meta_sub[, pc_cols_use] <- pc_scores_sub
    } else {
      pc_cols_use <- character(0)
    }
  } else {
    pc_cols_use <- character(0)
  }
  
  # Scale PCs AFTER they exist
  scaled_pc_names <- paste0("scaled_", pc_cols_use)
  for (i in seq_along(pc_cols_use)) {
    meta_sub[[ scaled_pc_names[i] ]] <- as.numeric(scale(meta_sub[[ pc_cols_use[i] ]]))
  }
  
  
  # ----- Convert categorical vars to factors -----
  for (v in colnames(meta_sub)) {
    if (is.character(meta_sub[[v]])) meta_sub[[v]] <- as.factor(meta_sub[[v]])
  }
  
  # ----- Filter samples & genes -----
  keep_samples <- complete.cases(meta_sub[, c(pc_cols_use, bio_covs), drop = FALSE])
  meta_sub <- meta_sub[keep_samples, , drop = FALSE]
  counts <- counts[, meta_sub$unique_id, drop = FALSE]
  
  # Compute qc metrics
  # qc_metrics[[ct]] <- data.frame(
  #   celltype = ct,
  #   mean_libsize = mean(colSums(counts)),
  #   dropout_rate = sum(counts == 0) / length(counts),
  #   median_cpm = median(cpm_edgeR(counts, log=TRUE), na.rm=TRUE)
  # )
  
  
  # ######## Remove counts associated to x and y sex chromosomes ##############
  # # Load Ensembl database
  # edb <- EnsDb.Hsapiens.v86
  # 
  # # Clean Ensembl IDs (remove version suffix if present)
  # ensembl_ids <- gsub("\\..*$", "", rownames(counts))
  # 
  # # Get chromosome annotation
  # gene_annot <- AnnotationDbi::select(
  #   edb,
  #   keys     = ensembl_ids,
  #   keytype  = "GENEID",
  #   columns  = c("SEQNAME")
  # )
  # 
  # # Keep autosomal genes only
  # autosomal_genes <- gene_annot$GENEID[
  #   gene_annot$SEQNAME %in% as.character(1:22)
  # ]
  # 
  # autosomal_genes <- unique(autosomal_genes)
  # 
  # # Now subset the counts matrix
  # keep <- ensembl_ids %in% autosomal_genes
  # 
  # counts_autosomal <- counts[keep, , drop = FALSE]
  # 
  # message("Genes before filtering: ", nrow(counts))
  # message("Genes after removing sex chromosomes: ", nrow(counts_autosomal))
  # 
  # ### Save sex chromosome genes separately
  # sex_chr_genes <- gene_annot$GENEID[
  #   gene_annot$SEQNAME %in% c("X", "Y")
  # ]
  # 
  # sex_chr_genes <- unique(sex_chr_genes)
  # 
  # write.table(
  #   sex_chr_genes,
  #   file = "sex_chromosome_genes_removed.txt",
  #   quote = FALSE,
  #   row.names = FALSE,
  #   col.names = FALSE
  # )
  
  ### Create DGElist
  dge <- DGEList(counts = counts)
  
  # ----- Fit full model with all covariates -----
  full_terms <- c(scaled_pc_names, bio_covs)
  # design_full <- model.matrix(as.formula(paste("~1 +", paste(full_terms, collapse = " + "))), data = meta_sub)
  # 1. Build design (do NOT filter columns)
  design_full <- model.matrix(
    ~ scaled_PC1 + scaled_PC2 + scaled_PC3 + scaled_PC4 +
      ancestry + sex + scaled_age + scaled_BMI,
    data = meta_sub
  )
  
  stopifnot("(Intercept)" %in% colnames(design_full))
  
  
  ## --- design_full created above, e.g.
  ## design_full <- model.matrix(as.formula(paste("~", paste(full_terms, collapse=" + "))), data = meta_sub)
  # keep_genes <- rowSums(cpm_edgeR(dge) > 1) >= 0.25 * ncol(counts_autosomal)  # expressed in >=50% donors --> too aggressive filter? Maybe could try 25%
  # keep_genes <- rowSums(dge$counts > 0) > 0  # remove just 0 counts, which should be enough as a filter
  keep_genes <- filterByExpr(dge, design_full) # best practice, removes genes with insufficient power fi that specific contrast and genes that are weakly expressed in one sex (model-aware filter)
  dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
  if (nrow(dge) < 10) {
    message("Too few genes after filtering for ", ct)
    next
  }
  
  # Compute counts of tested genes x celltype
  saveRDS(list(total_genes_before_filter = nrow(counts), 
               genes_after_filter = nrow(dge)), 
          file.path(data_dir, paste0(ct, "_gene_counts_info.rds"))
  )
  # gene_counts_info <- readRDS("/home/ivm/DGE_results_LM/HSC_MPP/data/HSC_MPP_gene_counts_info.rds")
  
  
  
  
  
  # Ensure design_full is a matrix
  if (!is.matrix(design_full)) design_full <- as.matrix(design_full)
  
  # Helper to get a column safely from matrix/data.frame
  get_col_safe <- function(mat, colname) {
    if (!colname %in% colnames(mat)) return(NULL)
    # keep drop = FALSE to preserve matrix shape
    mat[, colname, drop = FALSE]
  }
  
  # 1) Check scaled_BMI existence and variability
  if ("scaled_BMI" %in% colnames(design_full)) {
    bmi_col <- get_col_safe(design_full, "scaled_BMI")
    # convert to vector for unique() check
    bmi_vec <- as.vector(bmi_col)
    n_non_na_unique <- length(unique(bmi_vec[!is.na(bmi_vec)]))
    if (n_non_na_unique <= 1) {
      message(" - scaled_BMI has <=1 unique non-NA value (", n_non_na_unique, "). Removing from design.")
      # remove column
      keep_cols <- setdiff(colnames(design_full), "scaled_BMI")
      design_full <- design_full[, keep_cols, drop = FALSE]
    } else {
      message(" - scaled_BMI retained (", n_non_na_unique, " unique non-NA values).")
    }
  } else {
    message(" - scaled_BMI not present in design.")
  }
  
  # 2) Remove any other constant / all-NA columns (columns with <=1 unique non-NA value)
  col_keep <- apply(design_full, 2, function(x) !all(is.na(x)))
  
  # Always keep intercept
  col_keep[colnames(design_full) == "(Intercept)"] <- TRUE
  
  design_full <- design_full[, col_keep, drop = FALSE]
  
  stopifnot("(Intercept)" %in% colnames(design_full))
  
  
  # 3) Sanity: ensure design has fewer params than samples (residual df > 2)
  n_samples <- nrow(design_full)  # model.matrix rows == samples
  n_params  <- ncol(design_full)
  resid_df <- n_samples - n_params
  message(" - samples = ", n_samples, "; params = ", n_params, "; resid df = ", resid_df)
  if (resid_df < 3) {
    message(" Skipping this covariate/model for ", ct, 
            " due to insufficient residual df (", resid_df, ").")
    next   # comment out while debugging and running code line by line
  }
  
  
  
  # design_full is now pruned and safe to use in lmFit()/voom()
  
  
  
  
  
  # # ----------Drop predictors that have no variance within a celltype
  # # --> removes problematic terms like BMI only where it cannot be estimated rather than globally
  # # identify columns with no variation (only 1 unique non-NA value)
  # const_cols <- sapply(design_full, function(x) length(unique(x[!is.na(x)])) <= 1)
  # 
  # # filter columns, NOT rows
  # design_full <- design_full[, !const_cols, drop = FALSE]
  
  
  # # ---------- OR Model BMI only when available
  # if (all(is.na(design_full$scaled_BMI))) {
  #   design_full <- design_full[, !colnames(design_full) %in% "scaled_BMI"]
  # }
  
  # Run linear model
  dge <- calcNormFactors(dge)
  
  # vobj <- voom(dge, design_full, plot = FALSE) # Transforms counts to log2-CPM and estimates mean-variance relationship to compute precision weights
  
  # Mean-variance trend from voom (using base R)
  pdf(
    file.path(plots_dir, paste0(ct, "_voom_mean_variance.pdf")),
    width = 6,
    height = 5
  )
  
  vobj <- voom(dge, design_full, plot = TRUE) # Transforms counts to log2-CPM and estimates mean-variance relationship to compute precision weights
  
  mtext(
    paste0(ct, " - all genes, all samples, not covariate-specific"),
    side = 3,
    line = 0.2,
    cex = 0.85
  )
  
  dev.off()
  
  
  fit_full <- lmFit(vobj, design_full)
  fit_full <- eBayes(fit_full)
  
  # Test coefficients
  topTable(fit_full, coef = "sexF")
  topTable(fit_full, coef = "ancestryPakistani")
  topTable(fit_full, coef = "scaled_BMI")
  topTable(fit_full, coef = "scaled_age")
  
  saveRDS(fit_full, file.path(data_dir, paste0(ct, "_fit_full.rds")))
  
  # Compute variance explained per covariate
  ss_total <- rowSums((vobj$E - rowMeans(vobj$E))^2)
  for (term in colnames(design)) {
    ss_term <- rowSums((design[, term] * coef(fit_full)[, term])^2)
    variance_explained_list[[ct]][[term]] <- mean(ss_term / ss_total, na.rm=TRUE)
  }
  
  # ----- Test each bio covariate by dropping it -----
  
  # bio <- "sex"
  
  get_coef_name <- function(bio, coef_names) {
    if (bio %in% coef_names) return(bio)
    hits <- grep(paste0(bio), coef_names, value = TRUE)
    if (length(hits) != 1) {
      stop("Ambiguous or missing coefficient for ", bio)
    }
    hits
  }
  library(EnhancedVolcano)
  library(grid)
  library(ggplot2)
  
  for (bio in bio_covs) {
    # Find coefficients corresponding to that covariate
    # First check coefficients names
    colnames(fit_full$coefficients)
    print(bio)
    # Map bio_covs to coefficient names
    coef_names <- colnames(fit_full$coefficients)
    
    #coef_to_test <- get_coef_name(bio, colnames(fit_full$coefficients))
    coef_to_test <- colnames(fit_full$coefficients)[grep(bio,colnames(fit_full$coefficients))]
    
    # Generate topTable  
    tt <- topTable(fit_full, coef = coef_to_test, number = Inf, sort.by = "P") # or sort by adj.P.Val?
    tt_summary <- summary(tt$logFC)
    
    saveRDS(tt, file.path(data_dir, paste0(ct, "_", bio, "_topTable.rds")))
    
    ################################################# PLOTS #############################################
    ######################################################################################################
    
    # # Volcano plot BASE R
    # pdf(file.path(plots_dir, paste0(ct, "_DGE_", bio, "_volcano2.pdf")))
    # print(plot_volcano(tt, test_name = bio, ct = ct))
    # dev.off()
    
    ### Publication-ready volcano plot
    
    # Ensure Helvetica everywhere
    theme_set(theme_classic(base_family = "Helvetica", base_size = 12))
    
    # Defensive checks
    stopifnot(
      "logFC" %in% colnames(tt),
      "P.Value" %in% colnames(tt)
    )
    
    # Transform rownames into colname of tt for labelling
    tt$gene_symbol <- rownames(tt)
    tt$gene_symbol <- as.character(tt$gene_symbol)
    
    # Ensure logFC and pvalue are numeric
    tt$logFC    <- as.numeric(as.character(tt$logFC))
    tt$P.Value  <- as.numeric(as.character(tt$P.Value))
    
    
    
    # Remove NA rows
    tt <- tt[complete.cases(tt[, c("logFC", "P.Value")]), ]
    pdf(
      file.path(plots_dir, paste0(ct, "_volcano_", bio, ".pdf")),
      width = 6.5,
      height = 6
    )
    volcanoplot <- EnhancedVolcano(
      tt,
      lab = tt$gene_symbol,
      x = "logFC",
      y = "P.Value",
      pCutoff = 0.05,
      FCcutoff = 1,          # single abs value, not range i.e., 1 instead of c(-1, 1)
      pointSize = 2.0,       # bigger points
      labSize = 4.0,         # bigger labels
      title = paste0(ct, " - ", bio),
      subtitle = NULL,
      caption = NULL
    )
    print(volcanoplot)
    dev.off()
    
    # plot MA(mean-average)plot and voom mean-variance plot in one pdf
    pdf(
      file = file.path(
        plots_dir,
        paste0(ct, "_DGE_", bio, "_MAplot.pdf")
      ),
      width = 8,
      height = 10
    )
    
    ## 1. MA plot for sex
    #print(colnames(fit_full$coefficients)) # to check for sex coefficient name
    MAplot <- limma::plotMA(
      fit_full,
      coef = coef_to_test,
      main = paste(ct, " - ", bio, " - MA plot"),
      ylim = c(-4, 4)
    )
    print(MAplot)
    dev.off()
    
    message("Completed contrast for ", bio)
  }
  
  
  
  # ----- Step 2: Nested loop for diseases -----
  
  # ---- Detect all binary disease columns in phenotypes ----
  
  # Optionally skip first N metadata columns (e.g., skip 2:3)
  # Remove rows containing -1 in any phenotype column (excluding unique_id)
  phenotypes_sub <- phenotypes[, -c(1,2,3, 348), drop = FALSE]
  
  
  # Define disease traits columns
  id_col <- "unique_id"   # explicitly define ID column name, which is the only non-disease col
  disease_cols <- setdiff(colnames(phenotypes_sub), id_col)
  
  message("Detected ", length(disease_cols), " binary disease phenotypes.")
  print(disease_cols)
  
  # Filter phenotypes for a minimum of 50 cases
  # Assume:
  # - phenotypes_sub: data.frame with first column = sample ID, rest = diseases
  # - id_col: name of the sample ID column
  # - pheno_cols: all disease columns (excluding ID)
  
  min_cases <- 50
  
  # Step 1: check number of cases per disease
  valid_diseases <- sapply(disease_cols, function(d) {
    sum(phenotypes_sub[[d]] == 1, na.rm = TRUE) >= min_cases
  })
  
  # Step 2: keep only valid diseases
  filtered_diseases <- disease_cols[valid_diseases]
  
  message("Keeping ", length(filtered_diseases), " diseases out of ", length(disease_cols))
  
  # Keep only sample ID and filtered diseases
  phenotypes_sub <- phenotypes_sub[, c(id_col, filtered_diseases)]
  
  # Run disease loop
  # disease_col <- "ATOPIC_DERM"
  
  for (disease_col in filtered_diseases) {
    
    message("Running DGE for disease: ", disease_col)
    
    # merge phenotype column
    pheno_sub <- phenotypes[, c("unique_id", disease_col), drop = FALSE]
    merged_pheno <- merge(meta_sub, pheno_sub, by = "unique_id", all.x = FALSE, all.y = FALSE)
    # convert to character
    merged_pheno$unique_id <- as.character(merged_pheno$unique_id)
    
    # keep only Control (0) and Case (1)
    merged_pheno <- merged_pheno[merged_pheno[[disease_col]] %in% c(0,1), , drop = FALSE]
    
    if(nrow(merged_pheno) < 4){
      message("Too few samples for disease ", disease_col, "; skipping.")
      next
    }
    
    # subset counts
    # Remove NAs in $unique_ID
    merged_pheno <- merged_pheno[!is.na(merged_pheno$unique_id), ]
    
    # Intersect IDs
    shared_ids <- intersect(merged_pheno$unique_id, colnames(counts))
    
    if (length(shared_ids) < 3) {
      message("Skipping: too few shared IDs")
      next
    }
    
    # Merge and subset counts based on safe IDs
    merged_pheno <- merged_pheno[match(shared_ids, merged_pheno$unique_id), , drop = FALSE]
    merged_pheno$unique_id <- as.character(merged_pheno$unique_id)
    
    saveRDS(merged_pheno, file.path(data_dir, paste0(ct, "_merged_pheno.rds")))
    
    counts_d <- counts[, shared_ids, drop = FALSE]
    
    # Check what was filtered out and why
    setdiff(merged_pheno$unique_id, colnames(counts))
    
    length(merged_pheno$unique_id)
    length(unique(merged_pheno$unique_id)) # no duplicates
    length(colnames(counts)) #counts number vs phenotypes length is there is mismatch
    setdiff(merged_pheno$unique_id, colnames(counts))[1:20]
    
    sum(is.na(merged_pheno$unique_id))
    head(merged_pheno$unique_id[is.na(merged_pheno$unique_id)])
    
    # # Compute case/control counts per disease
    # shared <- intersect(colnames(counts_list[[ct]]), merged_pheno$unique_id)
    # ph <- merged_pheno[merged_pheno$unique_id %in% shared, ]
    # 
    # counts_table <- data.frame(
    #   disease = filtered_diseases,
    #   cases = sapply(filtered_diseases, function(d) sum(ph[[d]] == 1, na.rm=TRUE)),
    #   controls = sapply(filtered_diseases, function(d) sum(ph[[d]] == 0, na.rm=TRUE)),
    #   donors = sapply(filtered_diseases, function(d) sum(!is.na(ph[[d]])))
    # )
    # 
    # case_counts[[ct]] <- counts_table
    
    
    # make disease a factor: Control=0, Case=1
    merged_pheno[[disease_col]] <- factor(merged_pheno[[disease_col]],
                                          levels = c(0,1),
                                          labels = c("Control","Case"))
    
    ### Create DGElist
    dge_d <- DGEList(counts = counts_d)
    
    
    # 1. Build design (do NOT filter columns)
    tab <- table(merged_pheno[[disease_col]], useNA = "ifany")
    print(tab)
    
    if (sum(tab >= 2) < 2) {
      message("Skipping ", disease_col, " in ", ct, ": insufficient cases/controls")
      next
    }
    
    
    
    # design matrix: PCs + bio covariates + disease column
    model_terms <- c(scaled_pc_names, bio_covs, disease_col)
    model_terms <- model_terms[model_terms %in% colnames(merged_pheno)]
    design_d <- model.matrix(as.formula(paste("~", paste(model_terms, collapse = " + "))),
                             data = merged_pheno)
    
    # # Can't specify all 300 diseases in matrix
    # design_d <- model.matrix(
    #   ~ scaled_PC1 + scaled_PC2 + scaled_PC3 + scaled_PC4 +
    #     ancestry + sex + scaled_age + scaled_BMI,
    #   data = merged_pheno
    # )
    
    stopifnot("(Intercept)" %in% colnames(design_d))
    
    
    ## --- design_full created above, e.g.
    ## design_full <- model.matrix(as.formula(paste("~", paste(full_terms, collapse=" + "))), data = meta_sub)
    # keep_genes <- rowSums(cpm_edgeR(dge) > 1) >= 0.25 * ncol(counts_autosomal)  # expressed in >=50% donors --> too aggressive filter? Maybe could try 25%
    # keep_genes <- rowSums(dge$counts > 0) > 0  # remove just 0 counts, which should be enough as a filter
    keep_genes_d <- filterByExpr(dge_d, design_d) # best practice, removes genes with insufficient power fi that specific contrast and genes that are weakly expressed in one sex (model-aware filter)
    dge_d <- dge_d[keep_genes_d, , keep.lib.sizes = FALSE]
    if (nrow(dge_d) < 10) {
      message("Too few genes after filtering for ", ct)
      next
    }
    
    # Ensure design_full is a matrix
    if (!is.matrix(design_d)) design_d <- as.matrix(design_d)
    
    # Helper to get a column safely from matrix/data.frame
    get_col_safe <- function(mat, colname) {
      if (!colname %in% colnames(mat)) return(NULL)
      # keep drop = FALSE to preserve matrix shape
      mat[, colname, drop = FALSE]
    }
    
    #### Fit model ###
    
    # Run linear model
    dge_d <- calcNormFactors(dge_d)
    
    vobj_d <- voom(dge_d, design_d, plot = FALSE) # Transforms counts to log2-CPM and estimates mean-variance relationship to compute precision weights
    
    # # Mean-variance trend from voom (using base R)
    # pdf(
    #   file.path(plots_dir, paste0(ct, disease_col, "_voom_mean_variance.pdf")),
    #   width = 6,
    #   height = 5
    # )
    # 
    # vobj_d <- voom(dge_d, design_d, plot = TRUE) # Transforms counts to log2-CPM and estimates mean-variance relationship to compute precision weights
    # 
    # mtext(
    #   paste0(ct, disease_col, " - all genes, all samples, not covariate-specific"),
    #   side = 3,
    #   line = 0.2,
    #   cex = 0.85
    # )
    # 
    # dev.off()
    
    
    fit_d <- lmFit(vobj_d, design_d)
    fit_d <- eBayes(fit_d)
    
    # find coefficient
    
    coef_name <- grep(
      paste0("^", disease_col),
      colnames(fit_d$coefficients),
      value = TRUE
    )
    if (length(coef_name) == 0) {
      message("Skipping ", ct, "  ", disease_col, ": coef not estimable")
      next
    }
    
    coef_name <- coef_name[1]
    
    
    
    
    
    
    
    
    
    # disease_coef <- grep(paste0("^", disease_col), colnames(design_d), value = TRUE)[1]
    # 
    # 
    # coef_idx <- which(colnames(fit_d$coefficients) == disease_coef)
    tt_d <- topTable(fit_d, coef = coef_name, number = Inf, sort.by = "P")
    tt_d$Gene <- rownames(tt_d)
    
    # save results
    saveRDS(tt_d, file.path(data_dir, paste0(ct, "_DGE_", disease_col, "_top.rds")))
    
    
    
    ################################################# PLOTS #############################################
    ######################################################################################################
    
    # # Volcano plot
    # pdf(file.path(plots_dir, paste0(ct, "_DGE_", disease_col, "_volcano.pdf")))
    # print(plot_volcano(tt_d, test_name = disease_col, ct = ct))
    # dev.off()
    # 
    
    ### Publication-ready volcano plot
    library(EnhancedVolcano)
    library(grid)
    library(ggplot2)
    
    # Ensure Helvetica everywhere
    theme_set(theme_classic(base_family = "Helvetica", base_size = 12))
    
    # Defensive checks
    stopifnot(
      "logFC" %in% colnames(tt_d),
      "P.Value" %in% colnames(tt_d)
    )
    
    # Transform rownames into colname of tt_d for labelling
    tt_d$gene_symbol <- rownames(tt_d)
    tt_d$gene_symbol <- as.character(tt_d$gene_symbol)
    
    # Ensure logFC and pvalue are numeric
    tt_d$logFC    <- as.numeric(as.character(tt_d$logFC))
    tt_d$P.Value  <- as.numeric(as.character(tt_d$P.Value))
    
    
    
    # Remove NA rows
    tt_d <- tt_d[complete.cases(tt_d[, c("logFC", "P.Value")]), ]
    
    # Run diagnostics
    cat("\nCELLTYPE:", ct,
        "\nTRAIT:", bio,
        "\nCOEF USED:", coef_name,
        "\n")
    
    cat("coef exists:",
        coef_name %in% colnames(fit_d$coefficients), "\n")
    
    cat("non-NA coefficients:",
        sum(!is.na(fit_d$coefficients[, coef_name])), "\n")
    
    cat("non-NA p-values:",
        sum(!is.na(fit_d$p.value[, coef_idx])), "\n")
    
    cat("residual df:",
        fit_d$df.residual[1], "\n")
    
    
    pdf(
      file.path(plots_dir, paste0(ct, "_volcano_", disease_col, ".pdf")),
      width = 6.5,
      height = 6
    )
    
    volcanoplot_d <- EnhancedVolcano(
      tt_d,
      lab = tt_d$gene_symbol,
      x = "logFC",
      y = "P.Value",
      pCutoff = 0.05,
      FCcutoff = 1,          # single abs value, not range i.e., 1 instead of c(-1, 1)
      pointSize = 2.0,       # bigger points
      labSize = 4.0,         # bigger labels
      title = paste0(ct, " - ", disease_col),
      subtitle = NULL,
      caption = NULL
    )
    
    print(volcanoplot_d)
    
    dev.off()
    
    
    # # plot MA(mean-average)plot
    # pdf(
    #   file = file.path(
    #     plots_dir,
    #     paste0(ct, "_DGE_", disease_col, "_MAplot.pdf")
    #   ),
    #   width = 8,
    #   height = 10
    # )
    
    # ## 1. MA plot for sex
    # print(colnames(fit_full$coefficients)) # to check for sex coefficient name
    # 
    # limma::plotMA(
    #   fit_full,
    #   coef = "sexF",
    #   main = paste(ct, disease_col, "MA plot"),
    #   ylim = c(-4, 4)
    # )
    
    # dev.off()
    
    
  }
  
} # close celltype loop

# Bind rows of case/control counts per disease
case_counts_df <- bind_rows(case_counts, .id="celltype")

# Validation case/counts matrix
case_counts_df %>%
  group_by(disease) %>%
  summarise(total_cases = sum(cases)) %>%
  arrange(total_cases)

case_counts_df %>%
  group_by(disease) %>%
  summarise(unique_cases = n_distinct(cases))


# Bind rows for computed variance explained across cell types
variance_df <- do.call(rbind, lapply(names(variance_explained_list), function(ct) {
  data.frame(
    celltype = ct,
    term = names(variance_explained_list[[ct]]),
    var_exp = unlist(variance_explained_list[[ct]])
  )
}))

# Bind rows for qc metrics
qc_df <- bind_rows(qc_metrics)

# Save all outputs
write.csv(variance_df, "/home/ivm/DGE_results_LM_publication/summary_variance_explained.csv", row.names=FALSE)
# write.csv(genes_tested_df, "/home/ivm/DGE_results_LM/summary_gene_counts.csv", row.names=FALSE)
write.csv(case_counts_df, "/home/ivm/DGE_results_LM_publication/summary_case_counts.csv", row.names=FALSE)
write.csv(qc_df, "/home/ivm/DGE_results_LM_publication/summary_qc_metrics.csv", row.names=FALSE)









######################## CONVERT .RDS FILES INTO .CSV ##########################
library(tools)
library(fs)     # for path handling (comes with RStudio Server)
library(stringr)

# Input root directory
input_root <- "/home/ivm/DGE_results_LM_publication/"

# Output root directory
output_root <- "/home/ivm/DGE_results_LM_csv/"

# Create output root if not exists
if (!dir.exists(output_root)) dir.create(output_root, recursive = TRUE)

# Find all RDS files recursively
rds_files <- list.files(
  path = input_root,
  pattern = "\\.rds$",
  full.names = TRUE,
  recursive = TRUE
)

message("Found ", length(rds_files), " RDS files.")

for (f in rds_files) {
  message("\nProcessing: ", f)
  
  # Read file safely
  obj <- try(readRDS(f), silent = TRUE)
  
  if (inherits(obj, "try-error")) {
    message("  L Cannot read file, skipping.")
    next
  }
  
  if (!is.data.frame(obj) || nrow(obj) == 0) {
    message("    Not a non-empty data.frame, skipping.")
    next
  }
  
  # Construct relative path (everything after input_root)
  rel_path <- str_replace(f, paste0("^", input_root, "/"), "")
  
  # Change extension from .rds  .csv
  rel_path_csv <- sub("\\.rds$", ".csv", rel_path)
  
  # Build full output path
  out_file <- file.path(output_root, rel_path_csv)
  
  # Ensure output subdirectories exist
  out_dir <- dirname(out_file)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
    message("  =Á Created: ", out_dir)
  }
  
  # Save CSV
  write.csv(obj, out_file, row.names = TRUE)
  
  message("   Saved: ", out_file)
}

message("\n< All RDS  CSV conversions complete!")





############################### Gene annotation and heatmap ####################
library(pheatmap)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(stringr)

# Upload one file to see structure
B_naive <- read.csv2(file = "/home/ivm/DGE_results_LM_csv/B_naive/data/B_naive_ancestry_topTable.csv", header = TRUE, sep = ",")
head(B_naive)

# Directories
input_root <- "DGE_results_LM_csv"
output_root <- "DGE_results_LM_csv_annotated"

if (!dir.exists(output_root)) dir.create(output_root, recursive = TRUE)

# List all CSV files recursively
csv_files <- list.files(
  path = input_root,
  pattern = "\\.csv$",
  full.names = TRUE,
  recursive = TRUE
)

# Exclude summary stats files
exclude_files <- c(
  "summary_case_counts.csv",
  "summary_qc_metrics.csv",
  "summary_variance_explained.csv"
)

csv_files <- csv_files[!basename(csv_files) %in% exclude_files]

# Prepare gene annotation map from EnsDb
edb <- EnsDb.Hsapiens.v86
gene_map <- genes(edb, return.type = "DataFrame") %>%
  as.data.frame() %>%
  select(gene_id, gene_name) %>%
  distinct()

names(gene_map) <- c("ensembl_id", "gene_name")

# Process each CSV file
for (f in csv_files) {
  
  message("\nProcessing: ", f)
  
  df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  
  # Detect if first column has no name
  first_colname <- colnames(df)[1]
  
  if (is.na(first_colname) || first_colname == "" || first_colname == "X") {
    message("First column has no valid name  treating as Ensembl IDs")
    colnames(df)[1] <- "ensembl_id"
  } else {
    message("Unexpected format: first column already has a name (", first_colname, ")")
    message("Skipping file.")
    next
  }
  
  # Remove version numbers from Ensembl IDs (ENSG....xx ENSG....)
  df$ensembl_id <- str_replace(df$ensembl_id, "\\.\\d+$", "")
  
  # Annotate
  df_annot <- df %>%
    left_join(gene_map, by = "ensembl_id")
  
  # Prepare output filename & directory structure
  rel_path <- str_replace(f, paste0("^", input_root, "/"), "")
  out_path <- file.path(output_root, rel_path)
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  
  out_path <- sub("\\.csv$", "_annotated.csv", out_path)
  
  write.csv(df_annot, out_path, row.names = FALSE)
  
  message("  Saved annotated file: ", out_path)
}

message("All CSVs annotated successfully.")












### -------------------------- Plot Heatmap with proportion of DEGs -----------------------
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(viridis)
library(tibble)

#----------------------------------------------------------
# 1. Load all CSVs and compute DEG proportions
#----------------------------------------------------------
main_dir <- "DGE_results_LM_csv_annotated"

ct_dirs <- list.dirs(main_dir, recursive = FALSE, full.names = TRUE)

deg_list <- list()

for (ct_dir in ct_dirs) {
  celltype <- basename(ct_dir)
  data_dir <- file.path(ct_dir, "data")  # note: "Data" folder
  
  if (!dir.exists(data_dir)) next
  csvs <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(csvs) == 0) next
  
  message("Processing cell type: ", celltype)
  
  for (csv in csvs) {
    # Extract trait name:
    fname <- basename(csv)
    
    # Remove celltype prefix if present
    trait <- sub(paste0("^", celltype, "(_DGE)?_"), "", fname)
    # Remove _topTable, _top, _annotated, .csv
    trait <- sub("(_topTable|_top)?(_annotated)?\\.csv$", "", trait)
    
    # Read CSV
    df <- read.csv(csv, stringsAsFactors = FALSE)
    
    # Skip if no rows
    if (nrow(df) == 0) next
    
    # Compute total genes and significant genes
    total_genes <- nrow(df)
    sig_genes   <- sum(df$adj.P.Val < 0.05, na.rm = TRUE)
    prop_deg    <- sig_genes / total_genes
    
    deg_list[[length(deg_list) + 1]] <- data.frame(
      celltype = celltype,
      trait = trait,
      total_genes = total_genes,
      sig_genes = sig_genes,
      prop_deg = prop_deg,
      stringsAsFactors = FALSE
    )
  }
}

# Combine into a single data.frame
deg_df <- bind_rows(deg_list)

#----------------------------------------------------------
# 3. Order columns: first bio covs, then diseases
#----------------------------------------------------------
bio_covs <- c("scaled_age", "scaled_BMI", "sex", "ancestry")

deg_df <- deg_df %>%
  mutate(trait_type = ifelse(trait %in% bio_covs, "bio_cov", "disease")) %>%
  arrange(trait_type, trait)

write.csv(deg_df, file = "/home/ivm/DGE_results_LM_csv_annotated/degs_summary_table_DiseasexCelltype.csv", row.names = FALSE)

#----------------------------------------------------------
# 4. Heatmap (pheatmap)
#----------------------------------------------------------

# Convert deg_df to a matrix: rows = celltypes, cols = traits
# Assuming deg_df has columns: celltype, trait, prop_deg
# Pivot to wide format -for prop degs
deg_wide_prop <- deg_df %>%
  pivot_wider(
    id_cols = celltype,        # unique row identifiers
    names_from = trait,        # columns are traits
    values_from = prop_deg,
    values_fill = 0             # fill missing combinations with 0
  )

# Now check for duplicates
anyDuplicated(deg_wide_prop$celltype)   # should be 0

write.csv(deg_wide_prop, file = "/home/ivm/DGE_results_LM_csv_annotated/degs_wide_prop_DiseasexCelltype.csv", row.names = FALSE)

# Convert to matrix
deg_mat_prop <- deg_wide_prop %>%
  column_to_rownames("celltype") %>%
  as.matrix()

# Check
dim(deg_mat_prop)
head(deg_mat_prop)

# Remove empty celltype
deg_mat_prop <- deg_mat_prop[rownames(deg_mat_prop) != "HSC_MPP", ]

gap_pos <- sum(colnames(deg_mat_prop) %in% bio_covs)


# Save
pdf("/home/ivm/DGE_results_LM_csv_annotated/heatmap_deg_proportion_all.pdf", width=16, height=7)
pheatmap(
  mat = deg_mat_prop[,-(52:53)], # to remove the first 4 biocovs whose much higher values that were skewing all the palette
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = viridis(50, direction = -1),
  breaks = seq(0, max(deg_mat_prop[,-(1:4)], na.rm=TRUE), length.out=51),
  fontsize = 8,
  fontsize_row = 9,
  fontsize_col = 8,
  angle_col = 45,
  cellwidth = 14,
  cellheight = 12,
  gaps_col = gap_pos,
  main = "Proportion of DEGs per Trait per Cell Type",
  border_color = "grey90"
)
dev.off()

# Save
pdf("/home/ivm/DGE_results_LM_csv_annotated/heatmap_deg_proportion_diseasesOnly.pdf", width=17, height=7)
pheatmap(
  mat = deg_mat_prop[,-c(1:4, 52:53)], # to remove the first 4 biocovs that were skewing all the palette
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = viridis(50, direction = -1),
  #breaks = seq(0, max(deg_mat[,-(1:4)], na.rm=TRUE), length.out=51),
  fontsize = 8,
  fontsize_row = 9,
  fontsize_col = 8,
  angle_col = 45,
  cellwidth = 14,
  cellheight = 12,
  #gaps_col = gap_pos,
  main = "Proportion of DEGs per Trait per Cell Type",
  border_color = "grey90"
)
dev.off()

# Save
pdf("/home/ivm/DGE_results_LM_csv_annotated/heatmap_deg_proportion_biocovsOnly.pdf", width=14, height=7)
pheatmap(
  mat = deg_mat_prop[, 1:4], # to remove the first 4 biocovs that were skewing all the palette
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = viridis(50, direction = -1),
  #breaks = seq(0, max(deg_mat[,-(1:4)], na.rm=TRUE), length.out=51),
  fontsize = 8,
  fontsize_row = 9,
  fontsize_col = 8,
  angle_col = 45,
  cellwidth = 14,
  cellheight = 12,
  #gaps_col = gap_pos,
  main = "Proportion of DEGs per Trait per Cell Type",
  border_color = "grey90"
)
dev.off()


# Pivot to wide format -for absolute number of DEGs
deg_wide_abs <- deg_df %>%
  pivot_wider(
    id_cols = celltype,        # unique row identifiers
    names_from = trait,        # columns are traits
    values_from = sig_genes,
    values_fill = 0             # fill missing combinations with 0
  )

# Now check for duplicates
anyDuplicated(deg_wide_abs$celltype)   # should be 0

# Convert to matrix
deg_mat_abs <- deg_wide_abs %>%
  column_to_rownames("celltype") %>%
  as.matrix()

# Check
dim(deg_mat_abs)
head(deg_mat_abs)

# Remove empty celltype
deg_mat_abs <- deg_mat_abs[rownames(deg_mat_abs) != "HSC_MPP", ]



# Now plot absolute number of DEGs
# Save
pdf("/home/ivm/DGE_results_LM_csv_annotated/heatmap_deg_absNumber_all.pdf", width=17, height=7)
pheatmap(
  mat = deg_mat_abs[,-(52:53)], # to remove the first 4 biocovs whose much higher values that were skewing all the palette
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = viridis(50, direction = -1),
  breaks = seq(0, max(deg_mat_abs[,-(1:4)], na.rm=TRUE), length.out=51),
  fontsize = 8,
  fontsize_row = 9,
  fontsize_col = 8,
  angle_col = 45,
  cellwidth = 14,
  cellheight = 12,
  gaps_col = gap_pos,
  main = "Absolute number of DEGs per Trait per Cell Type",
  border_color = "grey90"
)
dev.off()

# Save
pdf("/home/ivm/DGE_results_LM_csv_annotated/heatmap_deg_absNumber_diseasesOnly.pdf", width=14, height=7)
pheatmap(
  mat = deg_mat_abs[,-c(1:4, 52:53)], # to remove the first 4 biocovs that were skewing all the palette
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = viridis(50, direction = -1),
  #breaks = seq(0, max(deg_mat[,-(1:4)], na.rm=TRUE), length.out=51),
  fontsize = 8,
  fontsize_row = 9,
  fontsize_col = 8,
  angle_col = 45,
  cellwidth = 14,
  cellheight = 12,
  #gaps_col = gap_pos,
  main = "Absolute number of DEGs per Trait per Cell Type",
  border_color = "grey90"
)
dev.off()

# Save
pdf("/home/ivm/DGE_results_LM_csv_annotated/heatmap_deg_absNumber_biocovsOnly.pdf", width=14, height=7)
pheatmap(
  mat = deg_mat_abs[, 1:4], # to remove the first 4 biocovs that were skewing all the palette
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = viridis(50, direction = -1),
  #breaks = seq(0, max(deg_mat[,-(1:4)], na.rm=TRUE), length.out=51),
  fontsize = 8,
  fontsize_row = 9,
  fontsize_col = 8,
  angle_col = 45,
  cellwidth = 14,
  cellheight = 12,
  #gaps_col = gap_pos,
  main = "Absolute number of DEGs per Trait per Cell Type",
  border_color = "grey90"
)
dev.off()










########################### Plot volcanos for selected traits and cell types ##########################
library(EnhancedVolcano)

print(filtered_diseases) # for info about available diseases analysed after filtering for minimum number of cases (n>=50)

# Select traits for volcano plots - based on Giuditta's analyses
selected_traits <- c("AB1_INTESTINAL_INFECTIONS", "D3_ANAEMIA_IRONDEF", "E4_HYTHY_AI_STRICT", "E4_OBESITY", "E4_HYPERCHOL", "E4_LIPOPROT", "H8_EXTERNAL", "H8_OTHEREAR", "I9_HYPTENS", 
                     "I9_HYPTENSESS", "I9_IHD", "J10_ASTHMA_EXMORE", "J10_PNEUMONIA", "K11_CHOLELITH", "HYPOTHYROIDISM", "MDD", "T2D", "sex", "ancestry", "scaled_age", "scaled_BMI"
) 


# Directory where annotated DEG CSVs are saved
data_root <- "/home/ivm/DGE_results_LM_csv_annotated/" # adjust to your root folder
celltypes_all <- list.dirs(data_root, recursive = FALSE, full.names = FALSE)

# # Select list of cell types of interest 
# celltypes_of_interest <- c("B_memory_IGHMhigh", "T_CD4_CTL", "T_CD4_EM", "T_CD4_Treg", "T_MAIT")
# celltypes <- intersect(celltypes_all, celltypes_of_interest)

# ct <- "T_CD4_CTL"
# trait <- "D3_ANAEMIANAS"

for(ct in celltypes_all){
  print(ct)
  
  ct_dir <- file.path(data_root, ct, "data")   # make sure folder is 'data' exactly
  if(!dir.exists(ct_dir)) {
    message("Skipping, folder does not exist: ", ct_dir)
    next
  }
  
  csv_files_all <- list.files(ct_dir, pattern = "_top.*\\.csv$", full.names = TRUE)
  if(length(csv_files_all) == 0) {
    message("No files found in ", ct_dir)
    next
  }
  
  # plot_dir <- file.path(data_root, ct, "plots_FDR0.05_LogFC1")
  plot_dir <- file.path(data_root, ct, "plots_FDR0.05_LogFC0.5")
  dir.create(plot_dir, showWarnings = FALSE)
  
  # Subset list of csvs based on traits of interest
  csv_files <- csv_files_all [
    sapply(csv_files_all, function(f)
      any(sapply(selected_traits, function(tr)
        grepl(tr, basename(f), ignore.case = TRUE)
      ))
    )
  ]
  
  # csv <- "/home/ivm/DGE_results_LM_csv_annotated/T_CD4_CTL/data/T_CD4_CTL_DGE_D3_ANAEMIANAS_top_annotated.csv"
  
  for(csv in csv_files){
    fname <- basename(csv)
    
    print(fname)
    # Define trait
    ## Case 1: disease traits (unchanged)
    if (grepl("_DGE_", fname)) {
      
      trait <- fname
      trait <- sub("^.*?_DGE_", "", trait)
      trait <- sub("_top.*$", "", trait)
      
    } else {
      
      ## Case 2: biological covariates
      ## extract trait by matching selected_traits
      trait <- selected_traits[
        sapply(selected_traits, function(t)
          grepl(paste0("_", t, "_"), fname))
      ]
      
      if (length(trait) == 0) {
        message("Skipping file (no matching trait): ", fname)
        next
      }
      
      trait <- trait[1]  # safety
    }
    
    
    # NOW trait = actual disease / bio-covariate name
    message("Processing: ", ct, " - ", trait)
    
    # OPTIONAL: filter to selected traits
    if (exists("selected_traits")) {
      if (!trait %in% selected_traits) next
    }
    
    
    
    tt <- read.csv(csv, stringsAsFactors = FALSE)
    
    # --- Force correct numeric types ---
    tt$logFC      <- as.numeric(tt$logFC)
    tt$P.Value  <- as.numeric(tt$P.Value)
    tt$adj.P.Val  <- as.numeric(tt$adj.P.Val)
    
    # Check if it worked or skip:
    if (!is.numeric(tt$logFC) | !is.numeric(tt$P.Value)) {
      message("Skipping (non-numeric): ", csv)
      next
    }
    
    
    # Optional: remove rows where conversion failed
    tt <- tt[!is.na(tt$logFC) & !is.na(tt$P.Value), ]
    
    # Skip if empty
    if (nrow(tt) == 0) {
      message("No rows in ", fname, " - skipping.")
      next
    }
    
    # Annotate significance column (p.value < 0.05)
    # tt <- tt %>%
    #   mutate(Sig = ifelse(adj.P.Val < 0.05, "Sig", "NS"))
    
    tt <- tt %>%
      mutate(
        SigAdj = adj.P.Val <= 0.05,
        SigRaw = P.Value < 0.05
      )
    
    
    if(!all(c("logFC", "P.Value", "adj.P.Val", "gene_name") %in% colnames(tt))) next
    
    
    
    # Print summary stats for every comparison  
    cat("\n====================================\n")
    cat("CELL TYPE:", ct, "\n")
    cat("TRAIT:", trait, "\n")
    cat("FILE:", csv, "\n")
    
    cat("N rows:", nrow(tt), "\n")
    
    cat("P.Value summary:\n")
    print(summary(tt$P.Value))
    
    cat("adj.P.Val summary:\n")
    print(summary(tt$adj.P.Val))
    
    cat("Sig p < 0.05:", sum(tt$P.Value < 0.05, na.rm=TRUE), "\n")
    cat("Sig FDR < 0.05:", sum(tt$adj.P.Val < 0.05, na.rm=TRUE), "\n")
    
    cat("logFC class:", class(tt$logFC), "\n")
    cat("P.Value class:", class(tt$P.Value), "\n")
    cat("adj.P.Val class:", class(tt$adj.P.Val), "\n")
    
    
    # Create argument for custom color dots
    keyvals <- ifelse(
      tt$SigAdj & abs(tt$logFC) >= 0.5, "Adj.P.Val & Log2FC", # Adjust if LogFC threshold changes
      ifelse(tt$SigAdj, "Adj.P.Val",
             ifelse(abs(tt$logFC) >= 0.5, "Log2FC", "NS")) # Adjust if LogFC threshold changes
    )
    
    keyvals <- factor(
      keyvals,
      levels = c("NS", "Log2FC", "Adj.P.Val", "Adj.P.Val & Log2FC")
    )
    
    cols <- c(
      "NS"        = "grey80",
      "Log2FC"   = "darkgreen",
      "Adj.P.Val" = "blue",
      "Adj.P.Val & Log2FC" = "red"
    )
    
    
    tt$gene_name <- as.character(tt$gene_name)
    
    top_genes <- tt %>% 
      dplyr::filter(!is.na(gene_name)) %>%
      dplyr::arrange(adj.P.Val) %>% 
      dplyr::slice_head(n = 15) %>% 
      dplyr::pull(gene_name)
    
    
    
    # Plot volcano
    #pdf(file.path(plot_dir, paste0("volcano_", trait, ".pdf")), width = 7, height = 6)
    enVol <- EnhancedVolcano(tt,
                             lab = tt$gene_name,
                             x = 'logFC',
                             y = 'P.Value',
                             pCutoff = 0.05,
                             FCcutoff = 0.5,
                             selectLab = top_genes,
                             title = paste0("DGE by ", trait, " - ", ct),
                             subtitle = "Raw p-values (y), FDR-based significance",
                             # caption = paste0("Cases: ", n_cases, " | Controls: ", n_ctrl),,
                             pointSize = 2.2,
                             labSize = 3.0,
                             colCustom = cols[keyvals],
                             drawConnectors = TRUE, 
                             widthConnectors = 0.5, 
                             colAlpha = 0.9,
                             legendPosition = 'right',
                             legendLabSize = 10,
                             legendIconSize = 3.0)
    ggsave(enVol,file=paste0(plot_dir, paste0("/volcano_", trait, ".pdf")),device="pdf",width = 8, height = 5)
    #dev.off()
  }
}



########################### Plot volcanos for bio covariates with high LFC threshold ##########################
library(EnhancedVolcano)

# Select traits for volcano plots
selected_traits <- c("sex", "ancestry", "scaled_BMI", "scaled_age" ) 


# Directory where annotated DEG CSVs are saved
data_root <- "/home/ivm/DGE_results_LM_csv_annotated/" # adjust to your root folder
celltypes_all <- list.dirs(data_root, recursive = FALSE, full.names = FALSE)

# ct <- "T_CD4_CTL"
# trait <- "D3_ANAEMIANAS"

for(ct in celltypes_all){
  print(ct)
  
  ct_dir <- file.path(data_root, ct, "data")   # make sure folder is 'data' exactly
  if(!dir.exists(ct_dir)) {
    message("Skipping, folder does not exist: ", ct_dir)
    next
  }
  
  csv_files_all <- list.files(ct_dir, pattern = "_top.*\\.csv$", full.names = TRUE)
  if(length(csv_files_all) == 0) {
    message("No files found in ", ct_dir)
    next
  }
  
  plot_dir <- file.path(data_root, ct, "Plots")
  dir.create(plot_dir, showWarnings = FALSE)
  
  # Subset list of csvs based on traits of interest
  csv_files <- csv_files_all [
    sapply(csv_files_all, function(f)
      any(sapply(selected_traits, function(tr)
        grepl(tr, basename(f), ignore.case = TRUE)
      ))
    )
  ]
  
  # csv <- "/home/ivm/DGE_results_LM_csv_annotated/T_CD4_CTL/data/T_CD4_CTL_DGE_D3_ANAEMIANAS_top_annotated.csv"
  
  # safe helper to escape regex metacharacters in strings
  escape_regex <- function(x) {
    gsub("([][{}()+*^$|\\\\.?])", "\\\\\\1", x)
  }
  
  for(csv in csv_files){
    fname <- basename(csv)
    
    print(fname)
    # Define trait
    trait <- fname
    # trait <- sub("^.*?_DGE_", "", trait)        # remove everything before "_DGE_"
    trait <- sub("_top.*$", "", trait)          # remove "_top", "_topTable", "_top_annotated"
    # 3) if a list of celltypes is provided, remove the leading celltype (longest first)
    if (!is.null(celltypes_all) && length(celltypes_all) > 0) {
      # ensure unique and sort by decreasing length so longer names match first
      cts <- unique(as.character(celltypes_all))
      cts <- cts[order(nchar(cts), decreasing = TRUE)]
      # escape regex metacharacters
      cts_esc <- escape_regex(cts)
      pattern <- paste0("^(", paste(cts_esc, collapse = "|"), ")_")  # match at start then underscore
      if (grepl(pattern, trait, ignore.case = TRUE, perl = TRUE)) {
        trait <- sub(pattern, "", trait, ignore.case = TRUE, perl = TRUE)
      }
    } else {
      # fallback: remove first block up to first underscore (only if no celltypes provided)
      trait <- sub("^[^_]+_", "", trait)
    }
    
    # NOW trait = actual disease / bio-covariate name
    message("Processing: ", ct, " - ", trait)
    
    # OPTIONAL: filter to selected traits
    if (exists("selected_traits")) {
      if (!trait %in% selected_traits) next
    }
    
    tt <- read.csv(csv, stringsAsFactors = FALSE)
    
    # --- Force correct numeric types ---
    tt$logFC      <- as.numeric(tt$logFC)
    tt$P.Value  <- as.numeric(tt$P.Value)
    tt$adj.P.Val  <- as.numeric(tt$adj.P.Val)
    
    # Check if it worked or skip:
    if (!is.numeric(tt$logFC) | !is.numeric(tt$adj.P.Val)) {
      message("Skipping (non-numeric): ", csv)
      next
    }
    
    
    # Optional: remove rows where conversion failed
    tt <- tt[!is.na(tt$logFC) & !is.na(tt$P.Value), ]
    
    # Skip if empty
    if (nrow(tt) == 0) {
      message("No rows in ", fname, " - skipping.")
      next
    }
    
    # Annotate significance column (p.value < 0.05)
    tt <- tt %>%
      mutate(Sig = ifelse(P.Value < 0.05, "Sig", "NS"))
    
    if(!all(c("logFC", "adj.P.Val", "gene_name") %in% colnames(tt))) next
    
    top_genes <- tt %>% arrange(adj.P.Val) %>% slice_head(n = 15) %>% pull(gene_name)
    
    
    # Print summary stats for every comparison  
    cat("\n====================================\n")
    cat("CELL TYPE:", ct, "\n")
    cat("TRAIT:", trait, "\n")
    cat("FILE:", csv, "\n")
    
    cat("N rows:", nrow(tt), "\n")
    
    cat("P.Value summary:\n")
    print(summary(tt$P.Value))
    
    cat("adj.P.Val summary:\n")
    print(summary(tt$adj.P.Val))
    
    cat("Sig p < 0.05:", sum(tt$P.Value < 0.05, na.rm=TRUE), "\n")
    cat("Sig FDR < 0.05:", sum(tt$adj.P.Val < 0.05, na.rm=TRUE), "\n")
    
    cat("logFC class:", class(tt$logFC), "\n")
    cat("P.Value class:", class(tt$P.Value), "\n")
    cat("adj.P.Val class:", class(tt$adj.P.Val), "\n")
    
    
    # Plot volcano
    #pdf(file.path(plot_dir, paste0("volcano_", trait, ".pdf")), width = 7, height = 6)
    enVol <- EnhancedVolcano(tt,
                             lab = tt$gene_name,
                             x = 'logFC',
                             y = 'adj.P.Val',
                             pCutoff = 0.05,
                             FCcutoff = 1.5,
                             selectLab = top_genes,
                             title = paste0("DGE by ", trait, " - ", ct),
                             subtitle = NULL,
                             # caption = paste0("Cases: ", n_cases, " | Controls: ", n_ctrl),,
                             pointSize = 2.2,
                             labSize = 3.0,
                             colAlpha = 0.9,
                             legendPosition = 'right',
                             legendLabSize = 10,
                             legendIconSize = 3.0)
    ggsave(enVol,file=paste0(plot_dir, paste0("/volcano_", trait, ".pdf")),device="pdf",width = 8, height = 5)
    #dev.off()
  }
}



#####################################################################################################################################################


#----------------------------------------------------------
# 5. Compute case/control ratios for diseases
# (Requires metadata: phenotype file with donors)
#----------------------------------------------------------

base_dir <- "/home/ivm/DGE_results_LM_csv_annotated/"
#filtered_counts_list   # named list of pseudobulk matrices (genes x samples)
# master_table_F3 # metadata with:
# - sample_id (matching colnames of counts)
# - disease columns coded as 0/1



library(dplyr)

case_counts_list <- list()

ct <- "B_naive"

for (ct in names(filtered_counts_list)) {
  message("Processing counts for cell type: ", ct)
  
  # Iteratively create results folders for specific outputs
  # Create folder for this cell type
  ct_dir <- file.path(base_dir, ct)
  metrics_dir  <- file.path(ct_dir, "metrics")
  
  dir.create(ct_dir,  showWarnings = FALSE)
  dir.create(metrics_dir,  showWarnings = FALSE)
  
  ### Number of PCs to include as technical covariates
  # number of genetic/technical PCs to use instead of pool/tranche
  pc_scores_df <- as.data.frame(pca_list[[ct]]$x) # extract PCA score for each cell type
  pc_scores_df$unique_id <- rownames(pc_scores_df) #create unique_id column from rownames
  rownames(pc_scores_df) <- NULL
  n_PC_use <- 4
  pc_scores_df_name <- "pc_scores_df"   # data.frame with columns: unique_id, PC1..PCn
  
  # Input data
  counts_list <- filtered_counts_list    # genes x samples
  
  bio_covs <- c("ancestry", "scaled_age", "sex", "scaled_BMI")
  pc_prefix <- "PC"      # columns in pc_scores_df
  
  out_dir <- "DGE_results_LM"
  dir.create(out_dir, showWarnings = FALSE)
  
  counts <- counts_list[[ct]]
  if (is.null(counts) || !is.matrix(counts) || ncol(counts) < 4) {
    message("Skipping ", ct, ": invalid counts or too few samples.")
    next
  }
  # Remove NAs from $pool_id column in mastertable
  meta_all <- master_table_F3 [!is.na(master_table_F3$pool_id), ]
  
  # Combine pool_id and donor_id into a single column with an underscore
  meta_all$unique_id <- paste(meta_all$pool_id, meta_all$donor_id, sep = "_")
  
  # Check the first few to confirm
  head(meta_all$unique_id)
  
  # Align metadata
  shared_ids <- intersect(colnames(counts), meta_all$unique_id)
  counts <- counts[, shared_ids, drop = FALSE]
  meta_sub <- meta_all[match(shared_ids, meta_all$unique_id), , drop = FALSE]
  stopifnot(all(colnames(counts) == meta_sub$unique_id))
  
  # Scale continuous variables so that they all have mean = 0 and sd = 1, to improve model stability and make effect sizes comparable
  # When covariates have different numeric scales, scaling makes all covariates comparable in scale, helping model algorithm converge faster and more reliably
  meta_sub$scaled_age <- as.numeric(scale(meta_sub$age_at_recruitment_first))
  meta_sub$scaled_BMI <- as.numeric(scale(meta_sub$BMI))
  
  
  # ----- Attach PCs as technical covariates -----
  # Attach PCs
  if (exists(pc_scores_df_name, envir = .GlobalEnv)) {
    pc_scores_df <- get(pc_scores_df_name, envir = .GlobalEnv)
    pc_cols <- grep(paste0("^", pc_prefix), colnames(pc_scores_df), value = TRUE)
    if (length(pc_cols) > 0) {
      pc_cols_use <- head(pc_cols, n_PC_use)
      pc_scores_sub <- pc_scores_df[match(meta_sub$unique_id, pc_scores_df$unique_id), pc_cols_use, drop = FALSE]
      meta_sub[, pc_cols_use] <- pc_scores_sub
    } else {
      pc_cols_use <- character(0)
    }
  } else {
    pc_cols_use <- character(0)
  }
  
  # Scale PCs AFTER they exist
  scaled_pc_names <- paste0("scaled_", pc_cols_use)
  for (i in seq_along(pc_cols_use)) {
    meta_sub[[ scaled_pc_names[i] ]] <- as.numeric(scale(meta_sub[[ pc_cols_use[i] ]]))
  }
  
  
  # ----- Convert categorical vars to factors -----
  for (v in colnames(meta_sub)) {
    if (is.character(meta_sub[[v]])) meta_sub[[v]] <- as.factor(meta_sub[[v]])
  }
  
  # ----- Filter samples & genes -----
  keep_samples <- complete.cases(meta_sub[, c(pc_cols_use, bio_covs), drop = FALSE])
  meta_sub <- meta_sub[keep_samples, , drop = FALSE]
  counts <- counts[, meta_sub$unique_id, drop = FALSE]
  
  
  # ----- Filter phenotypes ----------
  # Optionally skip first N metadata columns (e.g., skip 2:3)
  # Remove rows containing -1 in any phenotype column (excluding unique_id)
  phenotypes_sub <- phenotypes[, -c(1,2,3, 348), drop = FALSE]
  
  
  # Define disease traits columns
  id_col <- "unique_id"   # explicitly define ID column name, which is the only non-disease col
  disease_cols <- setdiff(colnames(phenotypes_sub), id_col)
  
  message("Detected ", length(disease_cols), " binary disease phenotypes.")
  print(disease_cols)
  
  # Filter phenotypes for a minimum of 50 cases
  # Assume:
  # - phenotypes_sub: data.frame with first column = sample ID, rest = diseases
  # - id_col: name of the sample ID column
  # - pheno_cols: all disease columns (excluding ID)
  
  min_cases <- 50
  
  # Step 1: check number of cases per disease
  valid_diseases <- sapply(disease_cols, function(d) {
    sum(phenotypes_sub[[d]] == 1, na.rm = TRUE) >= min_cases
  })
  
  # Step 2: keep only valid diseases
  filtered_diseases <- disease_cols[valid_diseases]
  
  message("Keeping ", length(filtered_diseases), " diseases out of ", length(disease_cols))
  
  # Keep only sample ID and filtered diseases
  phenotypes_sub <- phenotypes_sub[, c(id_col, filtered_diseases)]
  
  
  ct_samples <- colnames(counts)
  
  # merge phenotype column
  # pheno_sub <- phenotypes[, c("unique_id", disease_col), drop = FALSE]
  merged_pheno <- merge(meta_sub, phenotypes_sub, by = "unique_id", all.x = FALSE, all.y = FALSE)
  # convert to character
  merged_pheno$unique_id <- as.character(merged_pheno$unique_id)
  
  # keep only Control (0) and Case (1)
  # merged_pheno <- merged_pheno[merged_pheno[[disease_col]] %in% c(0,1), , drop = FALSE]
  # 
  # if(nrow(merged_pheno) < 4){
  #   message("Too few samples for disease ", disease_col, "; skipping.")
  #   next
  # }
  
  # subset counts
  # Remove NAs in $unique_ID
  merged_pheno <- merged_pheno[!is.na(merged_pheno$unique_id), ]
  
  # Intersect IDs
  shared_ids <- intersect(merged_pheno$unique_id, colnames(counts))
  
  if (length(shared_ids) < 3) {
    message("Skipping: too few shared IDs")
    next
  }
  
  # Merge and subset counts based on safe IDs
  merged_pheno <- merged_pheno[match(shared_ids, merged_pheno$unique_id), , drop = FALSE]
  merged_pheno$unique_id <- as.character(merged_pheno$unique_id)
  
  # saveRDS(merged_pheno, file.path(data_dir, paste0(ct, "_merged_pheno.rds")))
  
  counts_d <- counts[, shared_ids, drop = FALSE]
  
  # Check what was filtered out and why
  setdiff(merged_pheno$unique_id, colnames(counts))
  
  length(merged_pheno$unique_id)
  length(unique(merged_pheno$unique_id)) # no duplicates
  length(colnames(counts)) #counts number vs phenotypes lenght is there is mismatch
  setdiff(merged_pheno$unique_id, colnames(counts))[1:20]
  
  sum(is.na(merged_pheno$unique_id))
  head(merged_pheno$unique_id[is.na(merged_pheno$unique_id)])
  
  # Compute case/control counts per disease
  shared <- intersect(colnames(counts_list[[ct]]), merged_pheno$unique_id)
  ph <- merged_pheno[merged_pheno$unique_id %in% shared, ]
  
  counts_table <- data.frame(
    disease = filtered_diseases,
    cases = sapply(filtered_diseases, function(d) sum(ph[[d]] == 1, na.rm=TRUE)),
    controls = sapply(filtered_diseases, function(d) sum(ph[[d]] == 0, na.rm=TRUE)),
    donors = sapply(filtered_diseases, function(d) sum(!is.na(ph[[d]])))
  )
  
  # case_counts_df[[ct]] <- counts_table
  
  counts_table$celltype <- ct
  counts_table$case_control_ratio <- with(counts_table, cases / controls)
  
  case_counts_list[[ct]] <- counts_table
}

case_counts_df <- bind_rows(case_counts_list)

write.csv(
  case_counts_df,
  file = "case_control_counts_by_celltype_and_disease.csv",
  row.names = FALSE
)

# Convert case/ctrl ratio to log scale, for better scaling and interpretability
case_counts_df <- case_counts_df %>%
  mutate(
    case_control_ratio = cases / pmax(controls, 1),   # avoid /0
    log2_ratio = log2(case_control_ratio)
  )


# --------------- Plot dot plot per celltype ----------------
# pdf(file.path(base_dir,
#               "dotplot_log2_case_control_ratio_per_celltype.pdf"),
#     width = 16, height = 12)
# 
# ggplot(case_counts_df,
#        aes(x = disease,
#            y = log2_ratio)) +
#   geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
#   geom_point(size = 2.5, alpha = 0.8) +
#   facet_wrap(~ celltype, scales = "free_y") +
#   theme_bw(base_size = 11) +
#   theme(
#     axis.text.x  = element_text(angle = 45, hjust = 1, size = 5),
#     strip.text   = element_text(face = "bold"),
#     panel.grid.major.x = element_blank()
#   ) +
#   labs(
#     x = "Disease",
#     y = "log2(Case / Control)",
#     title = "Case/Control logRatio Faceted per Cell Type"
#   )
# 
# dev.off()




# ------------- Plot case/control log ratios per disease colored by celltype (non-faceted) -----------------
pdf(file.path(base_dir, "dotplot_log2_ratio_per_disease.pdf"),
    width = 16, height = 8)

ggplot(case_counts_df,
       aes(x = disease,
           y = log2_ratio,
           color = celltype)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_point(size = 2.8, alpha = 0.85) +
  scale_color_viridis_d() +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 11),
    legend.position = "right"
  ) +
  labs(
    x = "Disease",
    y = "log2(Case / Control)",
    color = "Cell type",
    title = "Case/Control logRatio per Disease - Coloured by Celltype"
  )

dev.off()

# ----------- Dotplot non-faceted per celltype
pdf(file.path(base_dir, "dotplot_log2_ratio_per_disease.pdf"),
    width = 16, height = 7)

ggplot(case_counts_df,
       aes(x = disease,
           y = log2_ratio,
           color = celltype)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_point(size = 2.8, alpha = 0.85) +
  scale_color_viridis_d() +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) + theme(
    strip.clip = "off",
    plot.margin = margin(10, 40, 50, 100)
  ) +
  labs(
    x = "Disease",
    y = "log2(Case / Control)",
    color = "Cell type",
    title = "Case/Control logRatio per Disease per Celltype"
  )

dev.off()

# Dotplot per celltype
# pdf(file.path(base_dir,
#               "dotplot_log2_case_control_ratio_faceted_by_celltype.pdf"),
#     width = 15, height = 10)
# 
# ggplot(case_counts_df,
#        aes(x = disease,
#            y = log2_ratio,
#            color = disease)) +
#   geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5) +
#   geom_point(size = 2.5, alpha = 0.85) +
#   facet_wrap(~ celltype, scales = "free_y") +
#   scale_color_viridis_d() +
#   theme_bw(base_size = 12) +
#   theme(
#     axis.text.x  = element_text(angle = 45, hjust = 1),
#     legend.position = "none",
#     strip.text = element_text(face = "bold")
#   ) +
#   labs(
#     x = "Disease",
#     y = "log2(Case / Control)",
#     title = "Case/Control Imbalance per Cell Type"
#   )
# 
# dev.off()


# Faceted dotplot colored by logRatio value
pdf(file.path(base_dir,
              "dotplot_log2_case_control_ratio_faceted_by_celltype.pdf"),
    width = 15, height = 10)
ggplot(case_counts_df,
       aes(x = disease,
           y = log2_ratio,
           color = abs(log2_ratio))) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 3) +
  scale_color_viridis_c() +
  facet_wrap(~ celltype) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) 


dev.off()


































