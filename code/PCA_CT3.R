
########## RE-RUN PCA AT CT3 GRANULARITY ###################
################## Assess covariate relevance for linear mixed model #######################
install.packages("edgeR")
install.packages("ggplot2")
install.packages("pheatmap")
install.packages("reshape2")

# Load libraries
library(edgeR)
library(ggplot2)
library(pheatmap)
library(reshape2)

# Load datafiles
# master table - donor level (each row is a cell)
master_table_F3 <- read.csv("/home/rstudio-server/F3_UKB_adata_obs_with_metadata.csv")      # Master table contains: age, sex, ethnicity, date at PBMC collection

# Filtered counts list
readRDS(filtered_counts_list, file = "/home/rstudio-server/filtered_counts_list.rds")

# ---- User parameters ----
nPC <- 15
min_cpm_keep <- 1
min_samples_keep <- 2
r2_threshold <- 0.05


# Specify covariates of interest
# skip_cols <- c("unique_id")
covariate_names <- c("sex", "smoking_status_combined", 
                     "bmi", "age"
                     # ,"pool_id", "tranche_id", "state",
                     # "age_at_recruitment_sq"
)

# covariate_names <- covariate_names[covariate_names %in% colnames(master_table_F3_donorL)]

covariate_names <- covariate_names[covariate_names %in% colnames(master_table_F3)]

### Load necesary files
# Load filtered counts list
filtered_counts_list <- readRDS("/home/ivm/CT3_DGE_analysis/filtered_counts_list.rds")
# master table - donor level (each row is a donor, already aggregated)
master_table_F3 <- read.table("/genesandhealth/red/DanielaZanotti/CARDINAL_db/v4/data/v4_master_table.tsv", header = TRUE, sep = "\t", check.names = FALSE)      # Master table contains: age, sex, ethnicity, date at PBMC collection
str(master_table_F3)
head(master_table_F3)


# ---- Helper functions ----
get_r2_for_covariate <- function(pc_scores_vec, covariate_vec){
  if(length(unique(na.omit(covariate_vec))) < 2) return(NA_real_)
  df <- data.frame(PC = pc_scores_vec, x = covariate_vec)
  m <- try(lm(PC ~ x, data = df), silent = TRUE)
  if(inherits(m, "try-error")) return(NA_real_)
  summary(m)$r.squared
}

get_p_for_covariate <- function(pc_scores_vec, covariate_vec){
  if(length(unique(na.omit(covariate_vec))) < 2) return(NA_real_)
  df <- data.frame(PC = pc_scores_vec, x = covariate_vec)
  m <- try(lm(PC ~ x, data = df), silent = TRUE)
  if(inherits(m, "try-error")) return(NA_real_)
  summary(m)$coefficients[2,4]  # p-value for covariate
}

# ---- Containers ----
pca_list <- list()
logCPM_list <- list()
r2_list <- list()
pval_list <- list()
pc_var_explained <- list()


# ct <- "B_exhausted"

# ---- Main loop over cell types ----
for(ct in names(filtered_counts_list)){
  message("Processing cell type: ", ct)
  counts <- filtered_counts_list[[ct]]
  
  # Skip empty or too-small cell types
  if(is.null(counts) || ncol(counts) < 2){
    message("  -> skipping, not enough samples")
    next
  }
  
  # Align metadata
  # Align metadata
  meta <- master_table_F3
  meta$donor_id_match <- gsub("--", "__", master_table_F3$donor_uid_tpd)
  head(meta$donor_id_match)
  
  meta <- meta[match(colnames(counts), meta$donor_id_match), ]
  missing_meta_idx <- which(is.na(meta$unique_id))
  if(length(missing_meta_idx) > 0){
    message("  -> dropping ", length(missing_meta_idx), " samples with missing metadata")
    keep_cols <- setdiff(colnames(counts), colnames(counts)[missing_meta_idx])
    counts <- counts[, keep_cols, drop = FALSE]
    meta <- master_table_F3[match(colnames(counts), master_table_F3$unique_id), ]
  }
  if(ncol(counts) < 2){
    message("  -> not enough samples after metadata filtering. Skipping ", ct)
    next
  }
  
  # ---- edgeR normalization ----
  dge <- DGEList(counts = counts)
  dge <- calcNormFactors(dge)
  
  keep_genes <- rowSums(edgeR::cpm(dge) > min_cpm_keep) >= min_samples_keep
  dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
  if(nrow(dge) < 10){
    message("  -> very few genes after filtering. Skipping ", ct)
    next
  }
  
  logCPM <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
  logCPM_list[[ct]] <- logCPM
  
  # ---- PCA ----
  pca <- prcomp(t(logCPM), scale. = TRUE, center = TRUE)
  pcs_to_use <- 1:min(nPC, ncol(pca$x))
  pca_list[[ct]] <- pca
  var_exp <- (pca$sdev^2)/sum(pca$sdev^2)
  pc_var_explained[[ct]] <- var_exp
  
  # Save pca_list
  saveRDS(pca_list, file = "/home/ivm/CT3_DGE_analysis/pca_list.rds")
  
  
  # ---- R² and p-value matrices ----
  r2_mat <- matrix(NA, nrow = length(covariate_names), ncol = length(pcs_to_use),
                   dimnames = list(covariate_names, paste0("PC", pcs_to_use)))
  pval_mat <- r2_mat
  
  for(i in seq_along(covariate_names)){
    covn <- covariate_names[i]
    covv <- meta[[covn]]
    for(j in seq_along(pcs_to_use)){
      pcj <- pcs_to_use[j]
      pc_scores_vec <- pca$x[, pcj]
      r2_mat[i,j] <- get_r2_for_covariate(pc_scores_vec, covv)
      pval_mat[i,j] <- get_p_for_covariate(pc_scores_vec, covv)
    }
  }
  
  r2_list[[ct]] <- r2_mat
  pval_list[[ct]] <- pval_mat
  
  # ---- R² heatmap with asterisks for significance ----
  sig_annot <- matrix("", nrow = nrow(r2_mat), ncol = ncol(r2_mat))
  sig_annot[pval_mat < 0.05] <- "*"
  sig_annot[pval_mat < 0.01] <- "**"
  sig_annot[pval_mat < 0.001] <- "***"
  
  # pdf(file.path(out_plots_dir, paste0(ct, "_PC_covariates_R2_heatmap.pdf")), width=8, height=6)
  # pheatmap(r2_mat,
  #          cluster_rows = FALSE,
  #          cluster_cols = FALSE,
  #          display_numbers = sig_annot,
  #          number_color = "black",
  #          color = colorRampPalette(c("white","steelblue"))(50),
  #          main = paste(ct, "PC ~ covariate R²"),
  #          fontsize_number = 10,
  #          angle_col = 0)
  # dev.off()
  
  # # ---- PCA scatter plots colored by covariates ----
  # scores <- as.data.frame(pca$x[, pcs_to_use, drop = FALSE])
  # scores$unique_id <- rownames(scores)
  # scores <- merge(scores, meta, by = "unique_id", sort = FALSE)
  # scores <- scores[match(rownames(pca$x), scores$unique_id), ]
  # 
  # pdf(file.path(out_plots_dir, paste0(ct, "_PC5_PC6_by_covariates.pdf")), width=6, height=5)
  # for(covn in covariate_names){
  #   covv <- scores[[covn]]
  #   if(is.numeric(covv)){
  #     p <- ggplot(scores, aes(PC5, PC6, color = covv)) +
  #       geom_point(size = 2) + labs(color = covn, title = paste(ct, "- colored by", covn)) +
  #       theme_minimal()
  #   } else {
  #     p <- ggplot(scores, aes(PC5, PC6, color = as.factor(covv))) +
  #       geom_point(size = 2) + labs(color = covn, title = paste(ct, "- colored by", covn)) +
  #       theme_minimal()
  #   }
  #   print(p)
  # }
  # dev.off()
  
  # ---- PCs exceeding R² threshold ----
  pcs_to_include <- which(colSums(r2_mat > r2_threshold, na.rm = TRUE) > 0)
  if(length(pcs_to_include) > 0){
    message("  -> PCs exceeding R² threshold: ", paste0("PC", pcs_to_include, collapse=", "))
  } else {
    message("  -> No PC exceeds R² threshold")
  }
  
  attr(r2_list[[ct]], "pcs_to_include") <- pcs_to_include
}

# ---- Summarize maximum R² across PCs for each covariate ----
max_r2_list <- list()
sig_labels_list <- list()

for (ct in names(r2_list)) {
  r2_mat <- r2_list[[ct]]
  pval_mat <- pval_list[[ct]]
  
  if (is.null(r2_mat) || all(is.na(r2_mat))) {
    warning("Skipping ", ct, ": R2 matrix is empty or NA.")
    next
  }
  
  # Replace non-finite values with NA
  r2_mat[!is.finite(r2_mat)] <- NA
  pval_mat[!is.finite(pval_mat)] <- NA
  
  # Compute max R² per covariate across PCs
  covariates_all <- rownames(r2_mat)
  max_r2 <- sapply(covariates_all, function(x) {
    vals <- r2_mat[x, ]
    if (all(is.na(vals))) return(0)
    max(vals, na.rm = TRUE)
  })
  
  # Compute best (lowest) p-value
  best_pval <- sapply(covariates_all, function(x) {
    vals <- pval_mat[x, ]
    if (all(is.na(vals))) return(1)
    min(vals, na.rm = TRUE)
  })
  
  # Significance stars (always present, but empty if not significant)
  sig_labels <- ifelse(best_pval < 0.001, "***",
                       ifelse(best_pval < 0.01, "**",
                              ifelse(best_pval < 0.05, "*", "")))
  
  # Store results
  max_r2_list[[ct]] <- data.frame(
    covariate = covariates_all,
    max_R2 = as.numeric(max_r2),
    significance = sig_labels,
    cell_type = ct
  )
}

# ---- Combine across cell types ----
combined_max_r2 <- do.call(rbind, max_r2_list)

# =9 Ensure all covariates appear for all cell types (fill missing with 0)
all_covariates <- unique(combined_max_r2$covariate)
all_celltypes <- unique(combined_max_r2$cell_type)

expanded <- expand.grid(covariate = all_covariates, cell_type = all_celltypes)
combined_max_r2 <- merge(expanded, combined_max_r2, all.x = TRUE)
combined_max_r2$max_R2[is.na(combined_max_r2$max_R2)] <- 0
combined_max_r2$significance[is.na(combined_max_r2$significance)] <- ""

# ---- Reshape for heatmap ----
heatmap_data <- reshape2::dcast(combined_max_r2, covariate ~ cell_type, value.var = "max_R2")
rownames(heatmap_data) <- heatmap_data$covariate
heatmap_data <- as.matrix(heatmap_data[,-1])

# Get significance stars matrix
sig_data <- reshape2::dcast(combined_max_r2, covariate ~ cell_type, value.var = "significance")
sig_data <- as.matrix(sig_data[,-1])
rownames(sig_data) <- rownames(heatmap_data)

# ---- Plot heatmap ----
pdf("/home/ivm/CT3_DGE_analysis/Max_R2_covariates_heatmap_allCovariates.pdf", width = 10, height = 8)
pheatmap(
  heatmap_data,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("white", "gold", "red"))(100),
  display_numbers = sig_data,
  number_color = "black",
  main = "Maximum R² per Covariate per Cell Type (all covariates shown)",
  fontsize_row = 10,
  fontsize_col = 10,
  angle_col = 45
)
dev.off()

message(" Saved heatmap: Max_R2_covariates_heatmap_allCovariates.pdf")




############### Create summary stats table for an example celltype to check for donor distribution for each covariate
# Choose one example cell type (or loop through all)
example_ct <- "T_CD8_naive"   # change to a cell type in your list

counts <- filtered_counts_list2[[example_ct]]
meta <- master_table_F3[match(colnames(counts), master_table_F3$unique_id), ]

# Keep only relevant covariates
covariate_names <- c("age_at_recruitment_first", "sex", "ancestry", "BMI", "state")
covariate_names <- covariate_names[covariate_names %in% colnames(meta)]
force_categorical <- c("sex")

# ---- Create a summary table ----
summary_list <- list()


for (covn in covariate_names) {
  covv <- meta[[covn]]
  
  # Force certain numeric covariates to be categorical
  if (covn %in% force_categorical) {
    covv <- as.factor(covv)
  }
  
  if (is.numeric(covv) && !covn %in% force_categorical) {
    df <- data.frame(
      covariate = covn,
      type = "numeric",
      category = NA,
      mean = mean(covv, na.rm = TRUE),
      sd = sd(covv, na.rm = TRUE),
      median = median(covv, na.rm = TRUE),
      min = min(covv, na.rm = TRUE),
      max = max(covv, na.rm = TRUE),
      count = sum(!is.na(covv)),
      percent = NA
    )
  } else {
    tbl <- as.data.frame(table(covv, useNA = "ifany"))
    colnames(tbl) <- c("category", "count")
    tbl$percent <- round(100 * tbl$count / sum(tbl$count), 1)
    df <- data.frame(
      covariate = covn,
      type = "categorical",
      category = tbl$category,
      mean = NA,
      sd = NA,
      median = NA,
      min = NA,
      max = NA,
      count = tbl$count,
      percent = tbl$percent
    )
  }
  summary_list[[covn]] <- df
}

# Combine results safely
common_cols <- Reduce(union, lapply(summary_list, names))
summary_list_aligned <- lapply(summary_list, function(x) {
  missing_cols <- setdiff(common_cols, names(x))
  for (mc in missing_cols) x[[mc]] <- NA
  x[common_cols]
})

#  Combine safely, even with mixed column structures
summary_table <- do.call(plyr::rbind.fill, summary_list_aligned)  # requires plyr package

# if you dont have plyr, use base R alternative:
# summary_table <- do.call(rbind, lapply(summary_list, function(x) {
#   common_cols <- Reduce(union, lapply(summary_list, names))
#   missing_cols <- setdiff(common_cols, names(x))
#   for (mc in missing_cols) x[[mc]] <- NA
#   x[common_cols]
# }))

# Save to file
out_file <- paste0("summary_covariates_", example_ct, ".tsv")
write.table(summary_table, out_file, sep = "\t", quote = FALSE, row.names = FALSE)

message(" Summary saved to: ", out_file)


#-----------------------------------------------------------------------------------------------------