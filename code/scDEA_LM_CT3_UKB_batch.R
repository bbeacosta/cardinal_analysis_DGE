#----------------------------------------- FREEZE 3 -----------------------------------------------
# scDEA_LM_CT3_UKB_batch.R 

##################################### LINEAR MODEL FOR DGE ###################################

# ---------------- PACKAGE SETUP (safer for batch jobs) ----------------
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(stringsAsFactors = FALSE)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Required packages (CRAN + Bioconductor)
pkgs <- c(
  "edgeR","limma","tidyr","dplyr","ggplot2","EnhancedVolcano",
  "optparse","stringi","tools","stringr","tibble","pheatmap",
  "EnsDb.Hsapiens.v86","AnnotationDbi","fs","viridis"
)

installed <- rownames(installed.packages())
missing <- setdiff(pkgs, installed)

if (length(missing) > 0) {
  message("Installing missing packages: ", paste(missing, collapse=", "))
  BiocManager::install(missing, ask = FALSE, update = FALSE)
}

# ggrepel pin only if missing or wrong version
if (!requireNamespace("ggrepel", quietly = TRUE) ||
    as.character(utils::packageVersion("ggrepel")) != "0.9.6") {
  remotes::install_version("ggrepel", version = "0.9.6", repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(edgeR); library(limma); library(tidyr); library(dplyr); library(ggplot2)
  library(EnhancedVolcano); library(optparse); library(stringr); library(stringi)
  library(tibble); library(pheatmap); library(EnsDb.Hsapiens.v86); library(AnnotationDbi)
  library(fs); library(viridis); library(tools)
})
# ---------------------------------------------------------------------

# ---------------- PATHS (single source of truth) ----------------
BASE_DIR <- "/home/rstudio-server"

MASTER_META_PATH <- file.path(BASE_DIR, "F3_UKB_adata_obs_with_metadata.csv")
PHENO_PATH       <- file.path(BASE_DIR, "disease_cases_controls_all.tsv")
FILTERED_PATH    <- file.path(BASE_DIR, "filtered_counts_list.rds")
PCA_PATH         <- file.path(BASE_DIR, "PCA_CT3", "pca_list.rds")

RESULTS_DGE_ROOT <- file.path(BASE_DIR, "results_DGE")
RESULTS_CSV_ROOT <- file.path(BASE_DIR, "results_DGE_csv")
RESULTS_ANN_ROOT <- file.path(BASE_DIR, "results_DGE_csv_annotated")

dir.create(BASE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RESULTS_DGE_ROOT, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(MASTER_META_PATH))
stopifnot(file.exists(PHENO_PATH))
stopifnot(file.exists(FILTERED_PATH))
stopifnot(file.exists(PCA_PATH))

# ---------------- LOAD INPUTS (read ONCE, correct formats) ----------------
master_table_F3 <- read.csv(MASTER_META_PATH, stringsAsFactors = FALSE, check.names = FALSE)
filtered_counts_list <- readRDS(FILTERED_PATH)
phenotypes <- read.delim(PHENO_PATH, stringsAsFactors = FALSE, check.names = FALSE)
pca_list <- readRDS(PCA_PATH)

# Keep only valid phenotype IDs
phenotypes <- phenotypes[!is.na(phenotypes$eid), ]

message("Loaded inputs:")
message(" - master_table_F3 rows: ", nrow(master_table_F3))
message(" - filtered_counts_list celltypes: ", length(filtered_counts_list))
message(" - phenotypes rows: ", nrow(phenotypes))
message(" - pca_list celltypes: ", length(pca_list))

# ================= HELPER FUNCTIONS =================
cpm_edgeR <- function(x, log = FALSE, prior.count = 0.25) {
  edgeR::cpm(x, log = log, prior.count = prior.count)
}

plot_volcano <- function(tt, test_name, ct, gene_label_col = "Gene", topN = 10) {
  df <- tt

  required_cols <- c("logFC", "P.Value", "adj.P.Val")
  if (!all(required_cols %in% colnames(df))) {
    stop("topTable object missing one of: logFC, P.Value, adj.P.Val")
  }

  if (!(gene_label_col %in% colnames(df))) {
    df[[gene_label_col]] <- rownames(df)
  }

  df$category <- dplyr::case_when(
    df$adj.P.Val < 0.05 & abs(df$logFC) >= 1 ~ "FDR < 0.05 & |LFC|",
    df$adj.P.Val < 0.05                     ~ "FDR < 0.05",
    abs(df$logFC) >= 1                      ~ "|LFC|",
    TRUE                                    ~ "Not significant"
  )

  colors <- c(
    "FDR < 0.05 & |LFC|" = "red",
    "FDR < 0.05"         = "blue",
    "|LFC|"              = "orange",
    "Not significant"    = "grey70"
  )

  top_labels <- df[order(df$adj.P.Val), ][1:min(topN, nrow(df)), ]

  ggplot(df, aes(x = logFC, y = -log10(P.Value), color = category)) +
    geom_point(alpha = 0.8, size = 1.2) +
    scale_color_manual(values = colors) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_text(
      data = top_labels,
      aes(label = !!rlang::sym(gene_label_col)),
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
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "right")
}

# variance/case/qc containers
variance_explained_list <- list()
case_counts <- list()
qc_metrics <- list()

# ---------------- IMPORTANT: ensure donor-level metadata exists ----------------
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


# ================= MAIN LOOP =================
for (ct in names(filtered_counts_list)) {
  message("\n===== CELL TYPE: ", ct, " =====")

  # Results dirs
  ct_dir    <- file.path(RESULTS_DGE_ROOT, ct)
  data_dir  <- file.path(ct_dir, "data")
  plots_dir <- file.path(ct_dir, "plots")
  dir.create(ct_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

  # Counts
  counts <- filtered_counts_list[[ct]]
  if (is.null(counts) || ncol(counts) < 4) {
    message("Skipping ", ct, ": invalid counts or too few samples.")
    next
  }

  # PCA scores for this ct (must exist)
  if (!ct %in% names(pca_list) || is.null(pca_list[[ct]]$x)) {
    message("Skipping ", ct, ": PCA scores missing in pca_list.")
    next
  }
  pc_scores_df <- as.data.frame(pca_list[[ct]]$x)
  pc_scores_df$unique_id <- rownames(pc_scores_df)
  rownames(pc_scores_df) <- NULL
  n_PC_use <- 4
  pc_prefix <- "PC"
  pc_cols <- grep(paste0("^", pc_prefix), colnames(pc_scores_df), value = TRUE)
  pc_cols_use <- head(pc_cols, n_PC_use)

  # Covariates expected in metadata
  bio_covs <- c("smoking_status_combined", "scaled_age", "sex", "scaled_bmi")

  # Metadata subset (remove NAs in pool_id if pool_id exists)
  meta_all <- master_table_F3_donorL_use
  if ("pool_id" %in% colnames(meta_all)) {
    meta_all <- meta_all[!is.na(meta_all$pool_id), , drop = FALSE]
  }

  # Required ID column to match counts colnames
  if (!("donor_uid_tpd_norm" %in% colnames(meta_all))) {
    stop("Metadata is missing required column donor_uid_tpd_norm (needed to align to counts colnames).")
  }

  # Align metadata to counts
  shared_ids <- intersect(colnames(counts), meta_all$donor_uid_tpd_norm)
  if (length(shared_ids) < 4) {
    message("Skipping ", ct, ": too few shared IDs between counts and metadata.")
    next
  }

  counts <- counts[, shared_ids, drop = FALSE]
  meta_sub <- meta_all[match(shared_ids, meta_all$donor_uid_tpd_norm), , drop = FALSE]
  stopifnot(all(colnames(counts) == meta_sub$donor_uid_tpd_norm))

  # Factors + scaling
  if ("sex" %in% colnames(meta_sub)) {
    meta_sub$sex <- factor(meta_sub$sex)
    if ("Male" %in% levels(meta_sub$sex)) meta_sub$sex <- relevel(meta_sub$sex, ref = "Male")
  } else {
    stop("Metadata missing 'sex' column.")
  }

  if ("smoking_status_combined" %in% colnames(meta_sub)) {
    meta_sub$smoking_status_combined <- factor(meta_sub$smoking_status_combined)
  } else {
    stop("Metadata missing 'smoking_status_combined' column.")
  }

  if (!("age" %in% colnames(meta_sub))) stop("Metadata missing 'age' column.")
  if (!("bmi" %in% colnames(meta_sub))) stop("Metadata missing 'bmi' column.")

  meta_sub$scaled_age <- as.numeric(scale(meta_sub$age))
  meta_sub$scaled_bmi <- as.numeric(scale(meta_sub$bmi))

  # Attach PCs (raw + scaled_PC*)
  if (length(pc_cols_use) > 0) {
    pc_scores_sub <- pc_scores_df[match(meta_sub$donor_uid_tpd_norm, pc_scores_df$unique_id), pc_cols_use, drop = FALSE]
    meta_sub[, pc_cols_use] <- pc_scores_sub
    for (pcn in pc_cols_use) {
      meta_sub[[paste0("scaled_", pcn)]] <- as.numeric(scale(meta_sub[[pcn]]))
    }
  } else {
    message("No PC columns detected for ", ct, " (expected PC1..).")
    next
  }

  # Convert character columns to factor
  for (v in colnames(meta_sub)) {
    if (is.character(meta_sub[[v]])) meta_sub[[v]] <- as.factor(meta_sub[[v]])
  }

  # Filter complete cases for terms used in the model
  scaled_pc_names <- paste0("scaled_", pc_cols_use)
  keep_samples <- complete.cases(meta_sub[, c(pc_cols_use, scaled_pc_names, bio_covs), drop = FALSE])
  meta_sub <- meta_sub[keep_samples, , drop = FALSE]
  counts   <- counts[, meta_sub$donor_uid_tpd_norm, drop = FALSE]

  if (ncol(counts) < 4) {
    message("Skipping ", ct, ": too few samples after filtering.")
    next
  }

  # Create DGEList
  dge <- DGEList(counts = counts)

  # Design matrix (CT3, 4 PCs)
  design_full <- model.matrix(
    ~ scaled_PC1 + scaled_PC2 + scaled_PC3 + scaled_PC4 +
      smoking_status_combined + sex + scaled_age + scaled_bmi,
    data = meta_sub
  )
  stopifnot("(Intercept)" %in% colnames(design_full))

  # Filter genes by expression for this design
  keep_genes <- filterByExpr(dge, design_full)
  dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
  if (nrow(dge) < 10) {
    message("Too few genes after filtering for ", ct)
    next
  }

  saveRDS(
    list(total_genes_before_filter = nrow(counts), genes_after_filter = nrow(dge)),
    file.path(data_dir, paste0(ct, "_gene_counts_info.rds"))
  )

  # Prune constant / all-NA design columns (keep intercept)
  col_ok <- apply(design_full, 2, function(x) !all(is.na(x)) && length(unique(x[!is.na(x)])) > 1)
  col_ok[colnames(design_full) == "(Intercept)"] <- TRUE
  design_full <- design_full[, col_ok, drop = FALSE]
  stopifnot("(Intercept)" %in% colnames(design_full))

  # Ensure residual df
  resid_df <- nrow(design_full) - ncol(design_full)
  message(" - samples = ", nrow(design_full), "; params = ", ncol(design_full), "; resid df = ", resid_df)
  if (resid_df < 3) {
    message("Skipping ", ct, ": insufficient residual df.")
    next
  }

  # Fit LM
  dge <- calcNormFactors(dge)

  pdf(file.path(plots_dir, paste0(ct, "_voom_mean_variance.pdf")), width = 6, height = 5)
  vobj <- voom(dge, design_full, plot = TRUE)
  mtext(paste0(ct, " - voom mean-variance"), side = 3, line = 0.2, cex = 0.85)
  dev.off()

  fit_full <- eBayes(lmFit(vobj, design_full))
  saveRDS(fit_full, file.path(data_dir, paste0(ct, "_fit_full.rds")))

  # ----- Test each bio covariate (t-test or F-test depending on factor levels) -----
  coef_names <- colnames(fit_full$coefficients)

  theme_set(theme_classic(base_family = "Helvetica", base_size = 12))

  for (bio in bio_covs) {
    message("Testing bio covariate: ", bio, " for ", ct)
    coef_to_test <- coef_names[grep(paste0("^", bio), coef_names)]
    if (length(coef_to_test) == 0) {
      message("No coefficients found for: ", bio, " -- skipping")
      next
    }

    if (length(coef_to_test) == 1) {
      tt <- topTable(fit_full, coef = coef_to_test, number = Inf, sort.by = "P")
      saveRDS(tt, file.path(data_dir, paste0(ct, "_", bio, "_topTable.rds")))
    } else {
      tt <- topTable(fit_full, coef = coef_to_test, number = Inf, sort.by = "F")
      saveRDS(tt, file.path(data_dir, paste0(ct, "_", bio, "_topTable_Ftest.rds")))
      for (one_coef in coef_to_test) {
        tt_one <- topTable(fit_full, coef = one_coef, number = Inf, sort.by = "P")
        saveRDS(tt_one, file.path(data_dir, paste0(ct, "_", one_coef, "_topTable.rds")))
      }
    }

    # Volcano plot only if it has the needed columns (and tt exists)
    if (exists("tt") && all(c("logFC", "P.Value") %in% colnames(tt))) {
      tt$gene_symbol <- rownames(tt)
      tt$logFC <- as.numeric(tt$logFC)
      tt$P.Value <- as.numeric(tt$P.Value)
      tt <- tt[complete.cases(tt[, c("logFC", "P.Value")]), , drop = FALSE]
      if (nrow(tt) > 0) {
        pdf(file.path(plots_dir, paste0(ct, "_volcano_", bio, ".pdf")), width = 6.5, height = 6)
        print(EnhancedVolcano(
          tt, lab = tt$gene_symbol, x = "logFC", y = "P.Value",
          pCutoff = 0.05, FCcutoff = 1,
          pointSize = 2.0, labSize = 4.0,
          title = paste0(ct, " - ", bio), subtitle = NULL, caption = NULL
        ))
        dev.off()
      }
    }

    # MA plot (use first coefficient if factor has multiple)
    ma_coef <- coef_to_test[1]
    pdf(file.path(plots_dir, paste0(ct, "_DGE_", bio, "_MAplot.pdf")), width = 8, height = 10)
    limma::plotMA(fit_full, coef = ma_coef, main = paste(ct, "-", bio, "- MA plot"), ylim = c(-4, 4))
    dev.off()

    message("Completed bio covariate: ", bio, " for ", ct)
  }

  # ----- Step 2: disease loop (binary traits with >= min_cases) -----
  phenotypes_sub <- phenotypes %>% dplyr::filter(!if_any(-eid, ~ . == -1))
  id_col <- "eid"
  disease_cols <- setdiff(colnames(phenotypes_sub), id_col)

  min_cases <- 50
  valid_diseases <- sapply(disease_cols, function(d) sum(phenotypes_sub[[d]] == 1, na.rm = TRUE) >= min_cases)
  filtered_diseases <- disease_cols[valid_diseases]
  message("Disease traits passing min_cases (", min_cases, ") for ", ct, ": ", length(filtered_diseases))

  if (length(filtered_diseases) == 0) next

  # Require eid to exist in metadata
  if (!("eid" %in% colnames(meta_sub))) {
    message("Skipping disease loop for ", ct, ": meta_sub has no 'eid' column.")
    next
  }

  for (disease_col in filtered_diseases) {
    message("Running DGE for disease: ", disease_col, " in ", ct)

    pheno_sub <- phenotypes[, c("eid", disease_col), drop = FALSE]
    merged_pheno <- merge(meta_sub, pheno_sub, by = "eid", all = FALSE)

    # keep Control/Case only
    merged_pheno <- merged_pheno[merged_pheno[[disease_col]] %in% c(0, 1), , drop = FALSE]
    if (nrow(merged_pheno) < 4) {
      message("Too few samples for disease ", disease_col, " in ", ct, "; skipping.")
      next
    }

    # subset counts
    shared_ids_d <- intersect(merged_pheno$donor_uid_tpd_norm, colnames(counts))
    if (length(shared_ids_d) < 4) {
      message("Skipping ", disease_col, " in ", ct, ": too few shared IDs.")
      next
    }

    merged_pheno <- merged_pheno[match(shared_ids_d, merged_pheno$donor_uid_tpd_norm), , drop = FALSE]
    counts_d <- counts[, shared_ids_d, drop = FALSE]

    merged_pheno[[disease_col]] <- factor(merged_pheno[[disease_col]],
                                          levels = c(0, 1),
                                          labels = c("Control", "Case"))

    tab <- table(merged_pheno[[disease_col]], useNA = "ifany")
    print(tab)
    if (sum(tab >= 2) < 2) {
      message("Skipping ", disease_col, " in ", ct, ": insufficient cases/controls")
      next
    }

    dge_d <- DGEList(counts = counts_d)

    # Build model terms
    model_terms <- c(scaled_pc_names, bio_covs, disease_col)
    model_terms <- model_terms[model_terms %in% colnames(merged_pheno)]
    design_d <- model.matrix(as.formula(paste("~", paste(model_terms, collapse = " + "))), data = merged_pheno)
    stopifnot("(Intercept)" %in% colnames(design_d))

    if (qr(design_d)$rank < ncol(design_d)) {
      message("Skipping ", ct, " - ", disease_col, ": design rank deficient")
      next
    }

    keep_genes_d <- filterByExpr(dge_d, design_d)
    dge_d <- dge_d[keep_genes_d, , keep.lib.sizes = FALSE]
    if (nrow(dge_d) < 10) {
      message("Too few genes after filtering for ", ct, " - ", disease_col)
      next
    }

    dge_d <- calcNormFactors(dge_d)
    vobj_d <- voom(dge_d, design_d, plot = FALSE)
    fit_d <- lmFit(vobj_d, design_d)

    if (all(fit_d$df.residual <= 0)) {
      message("Skipping ", ct, " - ", disease_col, ": no residual df")
      next
    }

    fit_d <- eBayes(fit_d)

    coef_name <- grep(paste0("^", disease_col), colnames(fit_d$coefficients), value = TRUE)
    if (length(coef_name) == 0) {
      message("Skipping ", ct, " - ", disease_col, ": coef not estimable")
      next
    }
    coef_name <- coef_name[1]

    tt_d <- topTable(fit_d, coef = coef_name, number = Inf, sort.by = "P")
    if (nrow(tt_d) == 0) {
      message("Skipping ", ct, " - ", disease_col, ": empty topTable")
      next
    }

    tt_d$Gene <- rownames(tt_d)
    saveRDS(tt_d, file.path(data_dir, paste0(ct, "_DGE_", disease_col, "_top.rds")))

    # Volcano
    if (all(c("logFC", "P.Value") %in% colnames(tt_d))) {
      tt_d$gene_symbol <- rownames(tt_d)
      tt_d$logFC <- as.numeric(tt_d$logFC)
      tt_d$P.Value <- as.numeric(tt_d$P.Value)
      tt_d <- tt_d[complete.cases(tt_d[, c("logFC", "P.Value")]), , drop = FALSE]
      if (nrow(tt_d) > 0) {
        pdf(file.path(plots_dir, paste0(ct, "_volcano_", disease_col, ".pdf")), width = 6.5, height = 6)
        print(EnhancedVolcano(
          tt_d,
          lab = tt_d$gene_symbol,
          x = "logFC",
          y = "P.Value",
          pCutoff = 0.05,
          FCcutoff = 1,
          pointSize = 2.0,
          labSize = 4.0,
          title = paste0(ct, " - ", disease_col),
          subtitle = NULL,
          caption = NULL
        ))
        dev.off()
      }
    }
  }
} # close celltype loop

message("\n=== DGE finished. Starting post-processing: RDS -> CSV, annotation, heatmaps, selected volcanos ===")

######################## CONVERT .RDS FILES INTO .CSV ##########################
input_root <- RESULTS_DGE_ROOT
output_root <- RESULTS_CSV_ROOT
if (!dir.exists(output_root)) dir.create(output_root, recursive = TRUE)

rds_files <- list.files(path = input_root, pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)
message("Found ", length(rds_files), " RDS files for CSV conversion.")

for (f in rds_files) {
  obj <- try(readRDS(f), silent = TRUE)
  if (inherits(obj, "try-error")) next
  if (!is.data.frame(obj) || nrow(obj) == 0) next

  rel_path <- stringr::str_replace(f, paste0("^", input_root, "/"), "")
  rel_path_csv <- sub("\\.rds$", ".csv", rel_path)
  out_file <- file.path(output_root, rel_path_csv)

  out_dir <- dirname(out_file)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  write.csv(obj, out_file, row.names = TRUE)
}

message("RDS -> CSV conversion complete.")

############################### Gene annotation ####################
input_root <- RESULTS_CSV_ROOT
output_root <- RESULTS_ANN_ROOT
if (!dir.exists(output_root)) dir.create(output_root, recursive = TRUE)

csv_files <- list.files(path = input_root, pattern = "\\.csv$", full.names = TRUE, recursive = TRUE)

exclude_files <- c("summary_case_counts.csv", "summary_qc_metrics.csv", "summary_variance_explained.csv")
csv_files <- csv_files[!basename(csv_files) %in% exclude_files]

edb <- EnsDb.Hsapiens.v86
gene_map <- genes(edb, return.type = "DataFrame") %>%
  as.data.frame() %>%
  dplyr::select(gene_id, gene_name) %>%
  dplyr::distinct()
names(gene_map) <- c("ensembl_id", "gene_name")

for (f in csv_files) {
  df <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)

  first_colname <- colnames(df)[1]
  if (is.na(first_colname) || first_colname == "" || first_colname == "X") {
    colnames(df)[1] <- "ensembl_id"
  } else {
    # already named -> likely not an Ensembl-id-first CSV; skip
    next
  }

  df$ensembl_id <- stringr::str_replace(df$ensembl_id, "\\.\\d+$", "")
  df_annot <- df %>% dplyr::left_join(gene_map, by = "ensembl_id")

  rel_path <- stringr::str_replace(f, paste0("^", input_root, "/"), "")
  out_path <- file.path(output_root, rel_path)
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  out_path <- sub("\\.csv$", "_annotated.csv", out_path)

  write.csv(df_annot, out_path, row.names = FALSE)
}

message("Annotation complete.")

### -------------------------- Heatmaps: DEG proportions + absolute counts -----------------------
main_dir <- RESULTS_ANN_ROOT
ct_dirs <- list.dirs(main_dir, recursive = FALSE, full.names = TRUE)

deg_list <- list()

for (ct_dir in ct_dirs) {
  celltype <- basename(ct_dir)
  data_dir <- file.path(ct_dir, "data")
  if (!dir.exists(data_dir)) next
  csvs <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(csvs) == 0) next

  for (csv in csvs) {
    fname <- basename(csv)
    trait <- sub(paste0("^", celltype, "(_DGE)?_"), "", fname)
    trait <- sub("(_topTable|_top)?(_annotated)?\\.csv$", "", trait)

    df <- read.csv(csv, stringsAsFactors = FALSE)
    if (nrow(df) == 0) next
    if (!("adj.P.Val" %in% colnames(df))) next

    total_genes <- nrow(df)
    sig_genes <- sum(df$adj.P.Val < 0.05, na.rm = TRUE)
    prop_deg <- sig_genes / total_genes

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

deg_df <- dplyr::bind_rows(deg_list)
if (nrow(deg_df) > 0) {
  bio_covs_hm <- c("scaled_age", "scaled_bmi", "sex", "smoking_status_combined")
  deg_df <- deg_df %>%
    dplyr::mutate(trait_type = ifelse(trait %in% bio_covs_hm, "bio_cov", "disease")) %>%
    dplyr::arrange(trait_type, trait)

  write.csv(deg_df, file = file.path(RESULTS_ANN_ROOT, "degs_summary_table_DiseasexCelltype.csv"), row.names = FALSE)

  deg_wide_prop <- deg_df %>%
    tidyr::pivot_wider(id_cols = celltype, names_from = trait, values_from = prop_deg, values_fill = 0)

  write.csv(deg_wide_prop, file = file.path(RESULTS_ANN_ROOT, "degs_wide_prop_DiseasexCelltype.csv"), row.names = FALSE)

  deg_mat_prop <- deg_wide_prop %>% tibble::column_to_rownames("celltype") %>% as.matrix()

  # Robust gaps position
  gap_pos <- sum(colnames(deg_mat_prop) %in% bio_covs_hm)

  pdf(file.path(RESULTS_ANN_ROOT, "heatmap_deg_proportion_all.pdf"), width = 40, height = 15)
  pheatmap::pheatmap(
    mat = deg_mat_prop,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = viridis::viridis(50, direction = -1),
    breaks = seq(0, max(deg_mat_prop, na.rm = TRUE), length.out = 51),
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

  deg_wide_abs <- deg_df %>%
    tidyr::pivot_wider(id_cols = celltype, names_from = trait, values_from = sig_genes, values_fill = 0)
  deg_mat_abs <- deg_wide_abs %>% tibble::column_to_rownames("celltype") %>% as.matrix()

  pdf(file.path(RESULTS_ANN_ROOT, "heatmap_deg_absNumber_all.pdf"), width = 40, height = 15)
  pheatmap::pheatmap(
    mat = deg_mat_abs,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = viridis::viridis(50, direction = -1),
    breaks = seq(0, max(deg_mat_abs, na.rm = TRUE), length.out = 51),
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
} else {
  message("No annotated DEG CSVs found for heatmaps (deg_df is empty).")
}

########################### Volcano plots for selected traits & cell types ##########################
selected_traits <- c(
  "AB1_INTESTINAL_INFECTIONS","D3_ANAEMIA_IRONDEF","E4_HYTHY_AI_STRICT","E4_OBESITY",
  "E4_HYPERCHOL","E4_LIPOPROT","H8_EXTERNAL","H8_OTHEREAR","I9_HYPTENS","I9_HYPTENSESS",
  "I9_IHD","J10_ASTHMA_EXMORE","J10_PNEUMONIA","K11_CHOLELITH","HYPOTHYROIDISM","MDD","T2D",
  "sex","ancestry","scaled_age","scaled_BMI"
)

data_root <- RESULTS_ANN_ROOT
celltypes_all <- list.dirs(data_root, recursive = FALSE, full.names = FALSE)

for (ct in celltypes_all) {
  ct_dir <- file.path(data_root, ct, "data")
  if (!dir.exists(ct_dir)) next

  csv_files_all <- list.files(ct_dir, pattern = "_top.*\\.csv$", full.names = TRUE)
  if (length(csv_files_all) == 0) next

  plot_dir <- file.path(data_root, ct, "plots_FDR0.05_LogFC1")
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

  csv_files <- csv_files_all[
    sapply(csv_files_all, function(f)
      any(sapply(selected_traits, function(tr) grepl(tr, basename(f), ignore.case = TRUE)))
    )
  ]

  for (csv in csv_files) {
    fname <- basename(csv)

    # trait parsing
    if (grepl("_DGE_", fname)) {
      trait <- sub("^.*?_DGE_", "", fname)
      trait <- sub("_top.*$", "", trait)
    } else {
      trait <- selected_traits[
        sapply(selected_traits, function(t) grepl(paste0("_", t, "_"), fname))
      ]
      if (length(trait) == 0) next
      trait <- trait[1]
    }

    if (!trait %in% selected_traits) next

    tt <- read.csv(csv, stringsAsFactors = FALSE)
    if (!all(c("logFC","P.Value","adj.P.Val","gene_name") %in% colnames(tt))) next

    tt$logFC <- as.numeric(tt$logFC)
    tt$P.Value <- as.numeric(tt$P.Value)
    tt$adj.P.Val <- as.numeric(tt$adj.P.Val)
    tt <- tt[!is.na(tt$logFC) & !is.na(tt$P.Value), , drop = FALSE]
    if (nrow(tt) == 0) next

    tt$gene_name <- as.character(tt$gene_name)

    keyvals <- ifelse(
      tt$adj.P.Val <= 0.05 & abs(tt$logFC) >= 1, "Adj.P.Val & Log2FC",
      ifelse(tt$adj.P.Val <= 0.05, "Adj.P.Val",
             ifelse(abs(tt$logFC) >= 1, "Log2FC", "NS"))
    )
    keyvals <- factor(keyvals, levels = c("NS","Log2FC","Adj.P.Val","Adj.P.Val & Log2FC"))

    cols <- c(
      "NS" = "grey80",
      "Log2FC" = "darkgreen",
      "Adj.P.Val" = "blue",
      "Adj.P.Val & Log2FC" = "red"
    )

    top_genes <- tt %>%
      dplyr::filter(!is.na(gene_name)) %>%
      dplyr::arrange(adj.P.Val) %>%
      dplyr::slice_head(n = 15) %>%
      dplyr::pull(gene_name)

    enVol <- EnhancedVolcano::EnhancedVolcano(
      tt,
      lab = tt$gene_name,
      x = "logFC",
      y = "P.Value",
      pCutoff = 0.05,
      FCcutoff = 1,
      selectLab = top_genes,
      title = paste0("DGE by ", trait, " - ", ct),
      subtitle = "Raw p-values (y), FDR-based significance",
      pointSize = 2.2,
      labSize = 3.0,
      colCustom = cols[keyvals],
      drawConnectors = TRUE,
      widthConnectors = 0.5,
      colAlpha = 0.9,
      legendPosition = "right",
      legendLabSize = 10,
      legendIconSize = 3.0
    )

    ggplot2::ggsave(
      plot = enVol,
      filename = file.path(plot_dir, paste0("volcano_", trait, ".pdf")),
      device = "pdf",
      width = 8,
      height = 5
    )
  }
}

message("\n=== DONE. Outputs written under:")
message(" - ", RESULTS_DGE_ROOT)
message(" - ", RESULTS_CSV_ROOT)
message(" - ", RESULTS_ANN_ROOT)
