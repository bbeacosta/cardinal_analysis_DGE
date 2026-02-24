
#----------------------------------------- FREEZE 3 -----------------------------------------------

##################################### LINEAR MIXED MODEL FOR DGE ###################################

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
# Load all necessary files for analyses
# load metadatafiles
# traits_macrocat <- read.csv("/genesandhealth/red/DanielaZanotti/data/disease/name_disease.csv", header = TRUE) # extended disease name for plotting and disease apparatus macrocategories

# master table - donor level (each row is a donor, already aggregated)
master_table_F3 <- read.table("/genesandhealth/red/DanielaZanotti/CARDINAL_db/v4/data/v4_master_table.tsv", header = TRUE, sep = "\t", check.names = FALSE)      # Master table contains: age, sex, ethnicity, date at PBMC collection
str(master_table_F3)
head(master_table_F3)

id_match_F3 <- read.table("/genesandhealth/red/DanielaZanotti/data/Freeze3/combined_status_v4_2.clean.tsv", header = TRUE, sep = "\t", check.names = FALSE)
str(id_match_F3)
head(id_match_F3)

# obs - cell level (each row is a cell, not aggregated by donor like pseudobulks)
obs_clean <- read.table("/genesandhealth/red/DanielaZanotti/data/Freeze3/obs_gh-ContaminationFree-noDoublets-annotated-QCed-counts_CLEAN.tsv", header = TRUE, sep = "\t", check.names = FALSE)   # Vacutainer ID
str(obs_clean)
head(obs_clean)

##### load Freeze 3 pseudobulk data #####
# Create counts directories
dir_bangpaki <- "/home/ivm/Desktop/Beatrice_Costa/gh-qced-cells/"

# extract everything before the first dot
all_files <- list.files(dir_bangpaki, full.names = FALSE)
celltypes <- unique(sub("\\..*", "", all_files))
print(celltypes) # 33 celltypes

# Set up ncell filter: at least 10 cells per donor per celltype (based on tables in: *ncells_x_donor.tsv)
# List all counts files
counts_files <- list.files(dir_bangpaki, pattern = "\\.raw\\.agg_sum\\.tsv$", full.names = TRUE)

# Output directory for filtered counts
out_dir <- file.path(dir_bangpaki, "filtered_counts")
dir.create(out_dir, showWarnings = FALSE)

# # Loop over counts files
# filtered_counts_list <- list()
# 
# for(cf in counts_files){
#   # ---- 1. Extract celltype name from filename
#   cf <- counts_files[6]
#   celltype <- str_replace(basename(cf), ".raw.agg_sum.tsv", "")
#   
#   message("Processing: ", celltype)
#   
#   # ---- 2. Match ncells file
#   ncells_file <- file.path(dir_bangpaki, paste0(celltype, ".ncells_x_donor.tsv"))
#   if(!file.exists(ncells_file)){
#     warning("No ncells file found for ", celltype, "  skipping")
#     next
#   }
#   
#   # 3. Read counts (genes x donors)
#   counts <- read.delim(cf, header=TRUE, sep="\t", stringsAsFactors=FALSE)
#   gene_col <- counts[[1]]
#   counts_mat <- counts[,-1, drop=FALSE]
#   rownames(counts_mat) <- gene_col
#   colnames(counts_mat) <- colnames(counts)[-1]
#   
#   # 4. Read ncells, filter rows with N_cells >= 10
#   ncells <- read.delim(ncells_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
#   ncells_keep <- ncells %>%
#     filter(N_cells >= 10)
#   
#   # 5. Keep donors corresponding to rows passing filter
#   donors_keep <- intersect(colnames(counts_mat), ncells_keep$unique.id)
#   
#   counts_filt <- counts_mat[, donors_keep, drop=FALSE]
#   message("Donors/columns kept: ", ncol(counts_filt), "/", ncol(counts_mat))
#   
#   
#   # 6. Save filtered counts
#   out_file <- file.path(out_dir, paste0(celltype, "_counts_filtered.tsv"))
#   write.table(
#     cbind(gene_id = rownames(counts_filt), counts_filt),
#     file=out_file,
#     sep="\t",
#     row.names=FALSE,
#     quote=FALSE
#   )
#   
#   # ---- 7. Store in list for R session use
#   filtered_counts_list[[celltype]] <- counts_filt
# }
# 
# # After loop: list of filtered counts per cell type
# names(filtered_counts_list)
# sapply(filtered_counts_list, ncol)
# 
# # Now drop celltypes with no donors >= 10 cells
# filtered_counts_list <- Filter(function(x) ncol(x) > 0, filtered_counts_list)
# sapply(filtered_counts_list, ncol)
# 
# # Save file for safe upload later
# saveRDS(filtered_counts_list, file = "/home/ivm/Desktop/Beatrice_Costa/results_DEA/filtered_counts_list.rds")

# ------------------------------------------------------------------------------

#Load counts 
filtered_counts_list <- readRDS("/home/ivm/Desktop/Beatrice_Costa/results_DEA/filtered_counts_list.rds")

# Create column in obs_clean that combines chromium_run_id and then _ donor_id columns to match pseudobulks
# Create a new column combining chromium_run_id and donor_id
obs_clean <- obs_clean %>%
  mutate(unique_id = paste(chromium_run_id, donor_id, sep = "_"))
colnames(obs_clean)
head(obs_clean$unique_id)
 
# Reduce obs_clean (cell-level metadata) to one row per donor (only keep unique ids for unique_id, state and vacutainer_id)
obs_id_map <- obs_clean %>%
  dplyr::select(vacutainer_id, unique_id, state) %>%
  distinct(vacutainer_id, .keep_all = TRUE)  # one per donor

# Join with master_table _F3 (donor-level metadata) --> vacutainer_id (obs_clean) with CARDINAL_ID_sample1 (mastertable) columns
master_table_F3 <- master_table_F3 %>%
  left_join(obs_id_map, by = c("CARDINAL_ID_sample1" = "vacutainer_id"))

# Keep only usable donors (follow $usable_donor_3, which Dani filtered for EU ancestry, frozen, duplicates)
# Keep only rows where usable_donor_3 is TRUE
master_table_F3 <- master_table_F3 %>%
  filter(usable_donor_3 == TRUE)

# Add age at collection squared to capture non-linear age effects
master_table_F3 <- master_table_F3 %>%
  mutate(age_at_recruitment_first_sq = age_at_recruitment_first^2)




######################### Work with pseudobulks ############################################

## --- Step 1: clean colnames in filtered_counts_list & build mapping ----
id_map_list <- lapply(names(filtered_counts_list), function(ct){
  mat <- filtered_counts_list[[ct]]
  old_ids <- colnames(mat)
  
  # extract substring starting at CRD_CM (remove string before that, so that matches with metadata)
  new_ids <- sub(".*(CRD_CM.*)", "\\1", old_ids)
  
  # Replace "__" with "_"
  new_ids <- gsub("__", "_", new_ids)
  
  # update colnames in the matrix
  colnames(mat) <- new_ids
  filtered_counts_list[[ct]] <<- mat   # update the global list
  
  # return mapping for all celltypes (elements in the counts list)
  data.frame(celltype   = ct,
             donor_id   = old_ids,
             unique_id  = new_ids,
             stringsAsFactors = FALSE)
})

# Combine mappings from all cell types into one dataframe
count_ids_df <- do.call(rbind, id_map_list)

## --- Step 2: merge with mastertable ----
id_match_adjusted <- merge(count_ids_df, master_table_F3, by = "unique_id")

# check overlap
cat("Number of matched unique_ids:", sum(count_ids_df$unique_id %in% master_table_F3$unique_id), "\n")

## --- Step 3: merge with bi_pheno ----
# merged_meta <- merge(id_match_adjusted, 
#                      bi_pheno, 
#                      by.x = "Vacutainer_ID", 
#                      by.y = "CARDINAL_ID_sample1")
# merged_full <- merge(merged_meta, covar,
#                      by.x = "Vacutainer_ID",
#                      by.y = "CARDINAL_ID_sample1")


# Flter count matrices by usable donors in master table
# Vector of usable donor IDs
usable_ids <- master_table_F3$unique_id

# Filter each count matrix in the list
filtered_counts_list2 <- lapply(filtered_counts_list, function(mat) {
  # Keep only columns (samples) present in usable donor list
  keep_cols <- intersect(colnames(mat), usable_ids)
  
  if (length(keep_cols) == 0) {
    warning("No usable donors found for a cell type; returning NULL.")
    return(NULL)
  }
  
  mat[, keep_cols, drop = FALSE]
})

# Optionally, remove empty cell types (if any became NULL)
filtered_counts_list2 <- Filter(Negate(is.null), filtered_counts_list2)

# Check results
sapply(filtered_counts_list2, ncol)   # number of donors per cell type after filtering

for(ct in names(filtered_counts_list2)){
  counts <- filtered_counts_list2[[ct]]
  cat("Cell type:", ct, "-> columns after filtering:", ncol(counts), "\n")
  
  # Confirm all colnames are in filtered_meta
  if(!all(colnames(counts) %in% master_table_F3$unique_id)){
    warning("Some counts columns not in filtered_meta for ", ct)
  }
}

# Confirm alignment for all celltypes in counts list with merged_full metadata
# Should return TRUE for all celltypes
for(ct in names(filtered_counts_list2)){
  counts <- filtered_counts_list2[[ct]]
  meta_ct <- master_table_F3[match(colnames(counts), master_table_F3$unique_id), ]
  
  # Check if all IDs match perfectly
  test <- all.equal(colnames(counts), meta_ct$unique_id)
  cat("Cell type:", ct, "-> alignment:", test, "\n")
}

##### Sanity check to check how many I filtered out and if numbers correspond
# Initialize results
state_summary_before <- list()
state_summary_after  <- list()

state_summary_before <- lapply(filtered_counts_list, function(mat) {
  ncol(mat)   # number of donors before filtering for this cell type
})

state_summary_after <- lapply(filtered_counts_list2, function(mat) {
  ncol(mat)   # number of donors after filtering for this cell type
})

# Combine and compare
before_df <- data.frame(
  celltype = names(state_summary_before),
  donors_before = unlist(state_summary_before)
)

after_df <- data.frame(
  celltype = names(state_summary_after),
  donors_after = unlist(state_summary_after)
)

summary_compare <- merge(before_df, after_df, by = "celltype", all = TRUE)
summary_compare$removed <- summary_compare$donors_before - summary_compare$donors_after

summary_compare




















# ################## Assess covariate relevance for linear mixed model #######################
# library(edgeR)
# library(ggplot2)
# library(pheatmap)
# library(reshape2)
# 
# # ---- User parameters ----
# nPC <- 15
# min_cpm_keep <- 1
# min_samples_keep <- 2
# r2_threshold <- 0.05
# out_plots_dir <- "PCA_plots"
# dir.create(out_plots_dir, showWarnings = FALSE)
# 
# # Specify covariates of interest
# skip_cols <- c("unique_id")
# covariate_names <- c("sex", "ancestry", "state",
#                      "age_at_recruitment", "BMI", 
#                      "pool_id", "tranche_id", 
#                      "age_at_recruitment_sq")
# 
# covariate_names <- covariate_names[covariate_names %in% colnames(master_table_F3)]
# 
# 
# # ---- Helper functions ----
# get_r2_for_covariate <- function(pc_scores_vec, covariate_vec){
#   if(length(unique(na.omit(covariate_vec))) < 2) return(NA_real_)
#   df <- data.frame(PC = pc_scores_vec, x = covariate_vec)
#   m <- try(lm(PC ~ x, data = df), silent = TRUE)
#   if(inherits(m, "try-error")) return(NA_real_)
#   summary(m)$r.squared
# }
# 
# get_p_for_covariate <- function(pc_scores_vec, covariate_vec){
#   if(length(unique(na.omit(covariate_vec))) < 2) return(NA_real_)
#   df <- data.frame(PC = pc_scores_vec, x = covariate_vec)
#   m <- try(lm(PC ~ x, data = df), silent = TRUE)
#   if(inherits(m, "try-error")) return(NA_real_)
#   summary(m)$coefficients[2,4]  # p-value for covariate
# }
# 
# # ---- Containers ----
# pca_list <- list()
# logCPM_list <- list()
# r2_list <- list()
# pval_list <- list()
# pc_var_explained <- list()
# 
# 
# 
# # ---- Main loop over cell types ----
# for(ct in names(filtered_counts_list2)){
#   message("Processing cell type: ", ct)
#   counts <- filtered_counts_list2[[ct]]
#   
#   # Skip empty or too-small cell types
#   if(is.null(counts) || ncol(counts) < 2){
#     message("  -> skipping, not enough samples")
#     next
#   }
#   
#   # Align metadata
#   meta <- master_table_F3[match(colnames(counts), master_table_F3$unique_id), ]
#   missing_meta_idx <- which(is.na(meta$unique_id))
#   if(length(missing_meta_idx) > 0){
#     message("  -> dropping ", length(missing_meta_idx), " samples with missing metadata")
#     keep_cols <- setdiff(colnames(counts), colnames(counts)[missing_meta_idx])
#     counts <- counts[, keep_cols, drop = FALSE]
#     meta <- master_table_F3[match(colnames(counts), master_table_F3$unique_id), ]
#   }
#   if(ncol(counts) < 2){
#     message("  -> not enough samples after metadata filtering. Skipping ", ct)
#     next
#   }
#   
#   # ---- edgeR normalization ----
#   dge <- DGEList(counts = counts)
#   dge <- calcNormFactors(dge)
#   
#   keep_genes <- rowSums(edgeR::cpm(dge) > min_cpm_keep) >= min_samples_keep
#   dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
#   if(nrow(dge) < 10){
#     message("  -> very few genes after filtering. Skipping ", ct)
#     next
#   }
#   
#   logCPM <- edgeR::cpm(dge, log = TRUE, prior.count = 1)
#   logCPM_list[[ct]] <- logCPM
#   
#   # ---- PCA ----
#   pca <- prcomp(t(logCPM), scale. = TRUE, center = TRUE)
#   pcs_to_use <- 1:min(nPC, ncol(pca$x))
#   pca_list[[ct]] <- pca
#   var_exp <- (pca$sdev^2)/sum(pca$sdev^2)
#   pc_var_explained[[ct]] <- var_exp
#   
#   # Save pca_list
#   saveRDS(pca_list, file = "pca_list.rds")
#   
#   
#   # ---- R² and p-value matrices ----
#   r2_mat <- matrix(NA, nrow = length(covariate_names), ncol = length(pcs_to_use),
#                    dimnames = list(covariate_names, paste0("PC", pcs_to_use)))
#   pval_mat <- r2_mat
#   
#   for(i in seq_along(covariate_names)){
#     covn <- covariate_names[i]
#     covv <- meta[[covn]]
#     for(j in seq_along(pcs_to_use)){
#       pcj <- pcs_to_use[j]
#       pc_scores_vec <- pca$x[, pcj]
#       r2_mat[i,j] <- get_r2_for_covariate(pc_scores_vec, covv)
#       pval_mat[i,j] <- get_p_for_covariate(pc_scores_vec, covv)
#     }
#   }
#   
#   r2_list[[ct]] <- r2_mat
#   pval_list[[ct]] <- pval_mat
#   
#   # ---- R² heatmap with asterisks for significance ----
#   sig_annot <- matrix("", nrow = nrow(r2_mat), ncol = ncol(r2_mat))
#   sig_annot[pval_mat < 0.05] <- "*"
#   sig_annot[pval_mat < 0.01] <- "**"
#   sig_annot[pval_mat < 0.001] <- "***"
#   
#   pdf(file.path(out_plots_dir, paste0(ct, "_PC_covariates_R2_heatmap.pdf")), width=8, height=6)
#   pheatmap(r2_mat,
#            cluster_rows = FALSE,
#            cluster_cols = FALSE,
#            display_numbers = sig_annot,
#            number_color = "black",
#            color = colorRampPalette(c("white","steelblue"))(50),
#            main = paste(ct, "PC ~ covariate R²"),
#            fontsize_number = 10,
#            angle_col = 0)
#   dev.off()
#   
#   # ---- PCA scatter plots colored by covariates ----
#   scores <- as.data.frame(pca$x[, pcs_to_use, drop = FALSE])
#   scores$unique_id <- rownames(scores)
#   scores <- merge(scores, meta, by = "unique_id", sort = FALSE)
#   scores <- scores[match(rownames(pca$x), scores$unique_id), ]
#   
#   pdf(file.path(out_plots_dir, paste0(ct, "_PC5_PC6_by_covariates.pdf")), width=6, height=5)
#   for(covn in covariate_names){
#     covv <- scores[[covn]]
#     if(is.numeric(covv)){
#       p <- ggplot(scores, aes(PC5, PC6, color = covv)) +
#         geom_point(size = 2) + labs(color = covn, title = paste(ct, "- colored by", covn)) +
#         theme_minimal()
#     } else {
#       p <- ggplot(scores, aes(PC5, PC6, color = as.factor(covv))) +
#         geom_point(size = 2) + labs(color = covn, title = paste(ct, "- colored by", covn)) +
#         theme_minimal()
#     }
#     print(p)
#   }
#   dev.off()
#   
#   # ---- PCs exceeding R² threshold ----
#   pcs_to_include <- which(colSums(r2_mat > r2_threshold, na.rm = TRUE) > 0)
#   if(length(pcs_to_include) > 0){
#     message("  -> PCs exceeding R² threshold: ", paste0("PC", pcs_to_include, collapse=", "))
#   } else {
#     message("  -> No PC exceeds R² threshold")
#   }
#   
#   attr(r2_list[[ct]], "pcs_to_include") <- pcs_to_include
# }
# 
# # ---- Summarize maximum R² across PCs for each covariate ----
# max_r2_list <- list()
# sig_labels_list <- list()
# 
# for (ct in names(r2_list)) {
#   r2_mat <- r2_list[[ct]]
#   pval_mat <- pval_list[[ct]]
#   
#   if (is.null(r2_mat) || all(is.na(r2_mat))) {
#     warning("Skipping ", ct, ": R2 matrix is empty or NA.")
#     next
#   }
#   
#   # Replace non-finite values with NA
#   r2_mat[!is.finite(r2_mat)] <- NA
#   pval_mat[!is.finite(pval_mat)] <- NA
#   
#   # Compute max R² per covariate across PCs
#   covariates_all <- rownames(r2_mat)
#   max_r2 <- sapply(covariates_all, function(x) {
#     vals <- r2_mat[x, ]
#     if (all(is.na(vals))) return(0)
#     max(vals, na.rm = TRUE)
#   })
#   
#   # Compute best (lowest) p-value
#   best_pval <- sapply(covariates_all, function(x) {
#     vals <- pval_mat[x, ]
#     if (all(is.na(vals))) return(1)
#     min(vals, na.rm = TRUE)
#   })
#   
#   # Significance stars (always present, but empty if not significant)
#   sig_labels <- ifelse(best_pval < 0.001, "***",
#                        ifelse(best_pval < 0.01, "**",
#                               ifelse(best_pval < 0.05, "*", "")))
#   
#   # Store results
#   max_r2_list[[ct]] <- data.frame(
#     covariate = covariates_all,
#     max_R2 = as.numeric(max_r2),
#     significance = sig_labels,
#     cell_type = ct
#   )
# }
# 
# # ---- Combine across cell types ----
# combined_max_r2 <- do.call(rbind, max_r2_list)
# 
# # =9 Ensure all covariates appear for all cell types (fill missing with 0)
# all_covariates <- unique(combined_max_r2$covariate)
# all_celltypes <- unique(combined_max_r2$cell_type)
# 
# expanded <- expand.grid(covariate = all_covariates, cell_type = all_celltypes)
# combined_max_r2 <- merge(expanded, combined_max_r2, all.x = TRUE)
# combined_max_r2$max_R2[is.na(combined_max_r2$max_R2)] <- 0
# combined_max_r2$significance[is.na(combined_max_r2$significance)] <- ""
# 
# # ---- Reshape for heatmap ----
# heatmap_data <- reshape2::dcast(combined_max_r2, covariate ~ cell_type, value.var = "max_R2")
# rownames(heatmap_data) <- heatmap_data$covariate
# heatmap_data <- as.matrix(heatmap_data[,-1])
# 
# # Get significance stars matrix
# sig_data <- reshape2::dcast(combined_max_r2, covariate ~ cell_type, value.var = "significance")
# sig_data <- as.matrix(sig_data[,-1])
# rownames(sig_data) <- rownames(heatmap_data)
# 
# # ---- Plot heatmap ----
# pdf("Max_R2_covariates_heatmap_allCovariates.pdf", width = 10, height = 8)
# pheatmap(
#   heatmap_data,
#   cluster_rows = TRUE,
#   cluster_cols = TRUE,
#   color = colorRampPalette(c("white", "gold", "red"))(100),
#   display_numbers = sig_data,
#   number_color = "black",
#   main = "Maximum R² per Covariate per Cell Type (all covariates shown)",
#   fontsize_row = 10,
#   fontsize_col = 10,
#   angle_col = 45
# )
# dev.off()
# 
# message(" Saved heatmap: Max_R2_covariates_heatmap_allCovariates.pdf")
# 
# 
# 
# 
# ############### Create summary stats table for an example celltype to check for donor distribution for each covariate
# # Choose one example cell type (or loop through all)
# example_ct <- "T_CD8_naive"   # change to a cell type in your list
# 
# counts <- filtered_counts_list2[[example_ct]]
# meta <- master_table_F3[match(colnames(counts), master_table_F3$unique_id), ]
# 
# # Keep only relevant covariates
# covariate_names <- c("age_at_recruitment_first", "sex", "ancestry", "BMI", "state")
# covariate_names <- covariate_names[covariate_names %in% colnames(meta)]
# force_categorical <- c("sex")
# 
# # ---- Create a summary table ----
# summary_list <- list()
# 
# 
# for (covn in covariate_names) {
#   covv <- meta[[covn]]
#   
#   # Force certain numeric covariates to be categorical
#   if (covn %in% force_categorical) {
#     covv <- as.factor(covv)
#   }
#   
#   if (is.numeric(covv) && !covn %in% force_categorical) {
#     df <- data.frame(
#       covariate = covn,
#       type = "numeric",
#       category = NA,
#       mean = mean(covv, na.rm = TRUE),
#       sd = sd(covv, na.rm = TRUE),
#       median = median(covv, na.rm = TRUE),
#       min = min(covv, na.rm = TRUE),
#       max = max(covv, na.rm = TRUE),
#       count = sum(!is.na(covv)),
#       percent = NA
#     )
#   } else {
#     tbl <- as.data.frame(table(covv, useNA = "ifany"))
#     colnames(tbl) <- c("category", "count")
#     tbl$percent <- round(100 * tbl$count / sum(tbl$count), 1)
#     df <- data.frame(
#       covariate = covn,
#       type = "categorical",
#       category = tbl$category,
#       mean = NA,
#       sd = NA,
#       median = NA,
#       min = NA,
#       max = NA,
#       count = tbl$count,
#       percent = tbl$percent
#     )
#   }
#   summary_list[[covn]] <- df
# }
# 
# # Combine results safely
# common_cols <- Reduce(union, lapply(summary_list, names))
# summary_list_aligned <- lapply(summary_list, function(x) {
#   missing_cols <- setdiff(common_cols, names(x))
#   for (mc in missing_cols) x[[mc]] <- NA
#   x[common_cols]
# })
# 
# #  Combine safely, even with mixed column structures
# summary_table <- do.call(plyr::rbind.fill, summary_list_aligned)  # requires plyr package
# 
# # if you dont have plyr, use base R alternative:
# # summary_table <- do.call(rbind, lapply(summary_list, function(x) {
# #   common_cols <- Reduce(union, lapply(summary_list, names))
# #   missing_cols <- setdiff(common_cols, names(x))
# #   for (mc in missing_cols) x[[mc]] <- NA
# #   x[common_cols]
# # }))
# 
# # Save to file
# out_file <- paste0("summary_covariates_", example_ct, ".tsv")
# write.table(summary_table, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
# 
# message(" Summary saved to: ", out_file)


#-----------------------------------------------------------------------------------------------------

#----------------------- RE-RUN LINEAR MIXED MODEL USING GENE EXPRESSION PCS ---------------------------------
# ====== PARAMETERS ======

# Convert all pseudobulk count tables to proper matrices
filtered_counts_list_fixed <- lapply(filtered_counts_list2, function(df) {
  # If it's a data frame with one column that contains nested lists or numeric vectors
  if (is.data.frame(df) && ncol(df) == 1 && is.list(df[[1]])) {
    mat <- do.call(cbind, df[[1]])
    colnames(mat) <- names(df[[1]])
    return(as.matrix(mat))
  }
  # If it's already a matrix
  if (is.matrix(df)) return(df)
  # If it's a data frame of numerics
  if (is.data.frame(df) && all(sapply(df, is.numeric))) {
    return(as.matrix(df))
  }
  return(NULL)
})

# Remove NULL or invalid entries
filtered_counts_list_fixed <- Filter(Negate(is.null), filtered_counts_list_fixed)

# Check structure
message(" After conversion, list has ", length(filtered_counts_list_fixed), " cell types.")
str(filtered_counts_list_fixed[[1]], max.level = 1)

saveRDS(filtered_counts_list_fixed, file = "/home/ivm/Desktop/Beatrice_Costa/results_DEA/filtered_counts_list_fixed.rds")



# Set wd
setwd("/home/ivm/LMM_DGE_results/")
outdir <- "/home/ivm/LMM_DGE_results/DGE_results_LRT_by_drop/"

#load PCA object
pca_list <- readRDS("/home/ivm/LMM_DGE_results/pca_list.rds")

# number of genetic/technical PCs to use instead of pool/tranche
pc_scores_df <- as.data.frame(pca_list[[ct]]$x) # extract PCA score for each cell type
pc_scores_df$unique_id <- rownames(pc_scores_df) #create unique_id column from rownames
rownames(pc_scores_df) <- NULL

n_PC_use <- 4
pc_scores_df_name <- "pc_scores_df"   # must be a data.frame with unique_id, PC1..PCn

# lists / data already in your environment
counts_list <- filtered_counts_list_fixed   # genes x samples (columns == unique_id)
meta_all <- master_table_F3                 # must contain unique_id, donor_id, and bio covariates

# biological covariates to test (tested by dropping one at a time)
bio_covs <- c("ancestry", "age_at_recruitment_first", "sex", "BMI")

# technical covariates replaced by PCs (do NOT include pool_id/tranche/state here)
pc_prefix <- "PC"   # assumed column names PC1, PC2, ...
random_effect <- "donor_id"

# output directory
out_dir <- "DGE_results_LRT_by_drop"
dir.create(out_dir, showWarnings = FALSE)

# parallel config for dream (adjust cores)
# register(SnowParam(workers = 4, type = "SOCK"))

# helper: explicit edgeR::cpm usage
cpm_edgeR <- function(x, log = FALSE, prior.count = 0.25) {
  edgeR::cpm(x, log = log, prior.count = prior.count)
}

# ====== MAIN LOOP ======

 ct <- "B_naive"
# for (ct in names(counts_list)) {
#   message("\n===== CELL TYPE: ", ct, " =====")
  counts <- counts_list[[ct]]
  if (is.null(counts) || !is.matrix(counts) || ncol(counts) < 4) {
    message("Skipping ", ct, ": invalid counts or too few samples.")
    next
  }
  
  # Align metadata by unique_id
  shared_ids <- intersect(colnames(counts), meta_all$unique_id)
  if (length(shared_ids) < 4) {
    message("Skipping ", ct, ": <4 overlapping samples between counts and metadata.")
    next
  }
  counts <- counts[, shared_ids, drop = FALSE]
  meta_sub <- meta_all[match(shared_ids, meta_all$unique_id), , drop = FALSE]
  stopifnot(all(colnames(counts) == meta_sub$unique_id))
  
  # drop unused factor levels
  meta_sub[] <- lapply(meta_sub, function(col){
    if (is.factor(col)) return(droplevels(col))
    if (is.character(col) && all(grepl("^[0-9.\\-]+$", na.omit(col)))) return(as.numeric(col))
    col
  })
  

  
  
  # ----- Attach PCs (technical) -----
  pcs_available <- exists(pc_scores_df_name, envir = .GlobalEnv)
  if (!pcs_available) {
    message("Warning: object '", pc_scores_df_name, "' not found. You must provide PCs in that data.frame to include PCs as technical covariates.")
    # Proceed without PCs, but recommended to supply them
    pc_terms <- character(0)
  } else {
    pc_scores_df <- get(pc_scores_df_name, envir = .GlobalEnv)
    # ensure unique_id present
    if (!"unique_id" %in% colnames(pc_scores_df)) {
      stop("pc_scores_df must contain column 'unique_id'")
    }
    # subset and order to meta_sub
    pcs_shared <- intersect(meta_sub$unique_id, pc_scores_df$unique_id)
    if (length(pcs_shared) < 1) {
      message("No PCs overlap with meta_sub; proceeding without PCs.")
      pc_terms <- character(0)
    } else {
      pc_scores_df2 <- pc_scores_df[match(meta_sub$unique_id, pc_scores_df$unique_id), , drop = FALSE]
      # check numeric PC columns exist
      pc_cols <- grep(paste0("^", pc_prefix), colnames(pc_scores_df2), value = TRUE)
      if (length(pc_cols) < 1) {
        message("No PC columns starting with '", pc_prefix, "' found in ", pc_scores_df_name, ". Proceeding without PCs.")
        pc_terms <- character(0)
      } else {
        # pick first n_PC_use PCs
        pc_cols_use <- head(pc_cols, n_PC_use)
        # add PCs into meta_sub as numeric columns
        meta_sub[, pc_cols_use] <- pc_scores_df2[, pc_cols_use, drop = FALSE]
        pc_terms <- pc_cols_use
        message("Added PCs: ", paste(pc_cols_use, collapse = ", "))
      }
    }
  }
  
  # ----- Select covariates to include in full model -----
  # Technical predictors = PCs (if present) else empty
  tech_terms <- pc_terms
  
  # Biological predictors (we will test them by dropping one at a time)
  bio_terms_present <- intersect(bio_covs, colnames(meta_sub))
  if (length(bio_terms_present) == 0) {
    message("No biological covariates found in metadata; skipping ", ct)
    next
  }
  
  # Convert char-like categorical columns to factor
  for (v in colnames(meta_sub)) {
    if (is.character(meta_sub[[v]])) meta_sub[[v]] <- as.factor(meta_sub[[v]])
  }
  
  # Drop any biological covariate that has <2 non-NA unique values in this subset (cannot test)
  bio_terms_present <- bio_terms_present[sapply(bio_terms_present, function(x) length(unique(na.omit(meta_sub[[x]]))) > 1)]
  if (length(bio_terms_present) == 0) {
    message("No variable biological covariates left for ", ct)
    next
  }
  
  # Combine full model terms (tech + all bio)
  full_terms <- c(tech_terms, bio_terms_present)
  # also ensure donor_id exists
  if (!random_effect %in% colnames(meta_sub)) {
    stop("random effect column '", random_effect, "' not found in metadata")
  }
  
  # Remove samples with NA in any of the predictors
  keep_samples <- complete.cases(meta_sub[, full_terms, drop = FALSE]) & !is.na(meta_sub[[random_effect]])
  if (sum(keep_samples) < 4) {
    message("Too few samples after dropping NAs (", sum(keep_samples), "). Skipping ", ct)
    next
  }
  meta_sub <- meta_sub[keep_samples, , drop = FALSE]
  counts <- counts[, meta_sub$unique_id, drop = FALSE]
  
  # Align metadata rows to counts columns
  # Keep only samples present in both counts and metadata
  common_ids <- intersect(colnames(counts), meta_sub$unique_id)
  
  # Subset counts
  counts <- counts[, common_ids, drop = FALSE]
  
  # Subset metadata and reorder to match counts columns
  meta_sub <- meta_sub[match(common_ids, meta_sub$unique_id), , drop = FALSE]
  
  # Set rownames of metadata to match counts columns
  rownames(meta_sub) <- meta_sub$unique_id
  
  # Verify alignment
  all(colnames(counts) == rownames(meta_sub))
  
  
  # ----- Collapse rare levels for factors (to avoid huge dummy expansion) -----
  # function to collapse rare levels (< min_count) into "Other"
  collapse_rare <- function(x, min_count = 3) {
    if (!is.factor(x)) return(x)
    tab <- table(x)
    rare <- names(tab)[tab < min_count]
    if (length(rare) == 0) return(droplevels(x))
    x2 <- as.character(x)
    x2[x2 %in% rare] <- "Other"
    return(as.factor(x2))
  }
  # apply to factor bio terms
  for (b in bio_terms_present) {
    if (is.factor(meta_sub[[b]])) meta_sub[[b]] <- collapse_rare(meta_sub[[b]], min_count = 3)
  }
  
  # ----- Sanity: Build model matrix and check rank -----
  full_formula <- as.formula(paste("~", paste(full_terms, collapse = " + "), "+ (1|", random_effect, ")"))
  
  # Drop factors with only one level (cannot use in contrasts)
  for (v in full_terms) {
    if (v %in% names(meta_sub) && is.factor(meta_sub[[v]])) {
      if (nlevels(droplevels(meta_sub[[v]])) < 2) {
        message(" Dropping ", v, " (only one level).")
        full_terms <- setdiff(full_terms, v)
      }
    }
  }
  
  
  # For building design matrix we'll ignore random effect
  mm <- model.matrix(as.formula(paste("~", paste(full_terms, collapse = " + "))), data = meta_sub, 
                     contrasts.arg = NULL)
  message("Model formula (fixed part): ", paste(full_terms, collapse = " + "))
  message("Design matrix dims: rows(samples)=", nrow(mm), " cols(parameters)=", ncol(mm))
  
  # If design has more columns than rows -> reduce technical PCs or drop rare factor levels further
  if (ncol(mm) >= nrow(mm)) {
    message("Design has >= samples. Attempting to reduce number of PCs to fit model.")
    # try reducing PCs gradually
    if (length(tech_terms) > 0) {
      for (k in seq(length(tech_terms)-1, 0, -1)) {
        try_terms <- c(head(tech_terms, k), bio_terms_present)
        mm_try <- model.matrix(as.formula(paste("~", paste(try_terms, collapse = " + "))), data = meta_sub)
        if (ncol(mm_try) < nrow(mm_try)) {
          message(" -> using ", k, " PCs instead of ", length(tech_terms))
          tech_terms <- head(tech_terms, k)
          full_terms <- c(tech_terms, bio_terms_present)
          mm <- mm_try
          break
        }
      }
    }
  }
  
  # If still rank-deficient, try removing highest-cardinality factor (common culprit: pool_id)
  if (ncol(mm) >= nrow(mm)) {
    # rank check
    if (qr(mm)$rank < ncol(mm)) {
      # find factor terms with most levels and drop them one by one
      factor_terms <- full_terms[sapply(full_terms, function(x) is.factor(meta_sub[[x]]))]
      if (length(factor_terms) > 0) {
        lvl_counts <- sapply(factor_terms, function(x) length(unique(meta_sub[[x]])))
        ord <- order(lvl_counts, decreasing = TRUE)
        for (i in ord) {
          dropt <- factor_terms[i]
          message("Dropping factor ", dropt, " (levels=", lvl_counts[i], ") to attempt identifiability")
          full_terms <- setdiff(full_terms, dropt)
          mm <- model.matrix(as.formula(paste("~", paste(full_terms, collapse = " + "))), data = meta_sub)
          if (ncol(mm) < nrow(mm) && qr(mm)$rank == ncol(mm)) break
        }
      }
    }
  }
  
  # final rank check
  if (ncol(mm) >= nrow(mm) || qr(mm)$rank < ncol(mm)) {
    message("Unable to form full-rank fixed model for ", ct, " (samples=", nrow(mm), ", params=", ncol(mm), "). Skipping.")
    next
  }
  
  # ----- Filter genes with explicit edgeR::cpm usage -----
  dge <- DGEList(counts = counts)
  keep_genes <- rowSums(edgeR::cpm(dge) > 1) >= 0.5 * ncol(dge) # test genes expressed in at least 50% of samples to make voom faster
  dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
  if (nrow(dge) < 10) {
    message("Too few genes after filtering for ", ct)
    next
  }
  dge <- calcNormFactors(dge)
  
  
  # ----- run voom + dream with automatic weights -----
  vobj <- voomWithDreamWeights(dge, full_formula, meta_sub, plot = FALSE)
  
  # ----- Fit full model with random effect using dream() -----
  fit_full <- tryCatch(
    dream(vobj, full_formula, meta_sub),
    error = function(e) {
      message("L dream failed for ", ct, ": ", e$message)
      return(NULL)
    }
  )
  if (is.null(fit_full)) next
  
  # Save full fit if wanted
  saveRDS(fit_full, file = file.path(out_dir, paste0(ct, "_fit_full.rds")))
  
  # ----- Test each biological covariate by dropping it and comparing via compareLRT -----
  for (bio in bio_terms_present) {
    message(" compareLRT test: drop ", bio)
    reduced_terms <- setdiff(full_terms, bio)
    if (length(reduced_terms) == 0) {
      message("  -> nothing left in reduced model; skipping")
      next
    }
    
    # ---- before building reduced_mm: ensure factors have no unused levels ----
    # (this should be done once after meta_sub is finalized for the celltype)
    facs <- sapply(meta_sub, is.factor)
    if (any(facs)) meta_sub[facs] <- lapply(meta_sub[facs], droplevels)
    
    # ensure numeric covariates are numeric (avoid accidental character columns)
    num_check <- sapply(meta_sub, is.numeric)
    # if some expected numeric columns are character, coerce carefully (example)
    # meta_sub$age_at_recruitment_first <- as.numeric(as.character(meta_sub$age_at_recruitment_first))
    
    reduced_mm <- model.matrix(as.formula(paste("~", paste(reduced_terms, collapse = " + "))), data = meta_sub)
    if (ncol(reduced_mm) >= nrow(reduced_mm) || qr(reduced_mm)$rank < ncol(reduced_mm)) {
      message("  -> reduced model not identifiable (params=", ncol(reduced_mm), " samples=", nrow(reduced_mm), "). Skipping LRT for ", bio)
      next
    }
    
    # double-check identifiability
    if (ncol(reduced_mm) >= nrow(reduced_mm) || qr(reduced_mm)$rank < ncol(reduced_mm)) {
      message("  -> reduced model not identifiable (params=", ncol(reduced_mm), " samples=", nrow(reduced_mm), "). Skipping LRT for ", bio)
      next
    }
    
    reduced_formula <- as.formula(paste("~", paste(reduced_terms, collapse = " + "), "+ (1|", random_effect, ")"))
    fit_reduced <- tryCatch({
      dream(vobj, reduced_formula, meta_sub)
    }, error = function(e) {
      message("  -> dream reduced failed for ", bio, ": ", e$message)
      return(NULL)
    })
    
    if (is.null(fit_reduced)) {
      message("L Reduced model failed for ", bio)
    } else {
      message(" Reduced model ok for ", bio)
    }
    

    if (is.null(fit_reduced)) next
  }

    # LRT (full vs reduced)
    # ----- Likelihood Ratio Test (manual) -----
    lrt <- tryCatch({
      # extract log-likelihoods and degrees of freedom
      ll_full <- sum(fit_full$logLik, na.rm = TRUE)
      ll_reduced <- sum(fit_reduced$logLik, na.rm = TRUE)
      df_full <- ncol(coef(fit_full))
      df_reduced <- ncol(coef(fit_reduced))
      
      
      # LRT statistic
      LRT_stat <- 2 * (ll_full - ll_reduced)
      df_diff <- df_full - df_reduced
      p_val <- pchisq(LRT_stat, df = df_diff, lower.tail = FALSE)
      
      data.frame(
        covariate = bio,
        LRT_stat = LRT_stat,
        df_diff = df_diff,
        p.value = p_val,
        FDR = p.adjust(p_val, method = "fdr")
      )
    }, error = function(e) {
      message("  -> LRT computation failed for ", bio, ": ", e$message)
      return(NULL)
    })
    
    # Extract summary table
    if (!is.null(lrt)) {
      saveRDS(list(
        full = fit_full,
        reduced = fit_reduced,
        LRT_summary = lrt
      ), file = file.path(out_dir, paste0(ct, "_drop_", bio, "_LRT.rds")))
      
      message("  -> LRT done for ", bio, ". p=", signif(lrt$p.value, 3), 
              " FDR=", signif(lrt$FDR, 3))
    }
  
  # Save bobj and meta_sub objects for future upload  
  saveRDS(vobj, file = file.path(outdir, paste0(ct, "_vobj.rds")))
  saveRDS(meta_sub, file = file.path(outdir, paste0(ct, "_meta_sub.rds")))
  
  
  #------------------------ MAKE PLOTS --------------------------------------
  ### Run your main dream model with fixed biological and random technical covariates (for DGE).
  ### Then separately run variancePartition with all covariates as random to visualize variance contributions.
  ### This is standard practice in multi-omic pseudobulk pipelines.
  outdir <- "/home/ivm/LMM_DGE_results/DGE_results_LRT_by_drop/"
  ct <- "B_naive"
  
  # Reload objects
  vobj <- readRDS(file.path(outdir, paste0(ct, "_vobj.rds")))
  meta_sub <- readRDS(file.path(outdir, paste0(ct, "_meta_sub.rds")))
  
  
  # ---- reuse vobj and meta_sub for this cell type ----
  # trasform covariates as factors in order to not have to force them to be random effects
  form_random <- ~ (1|ancestry) + sex + age_at_recruitment_first+ 
    BMI + PC1 + PC2 + PC3 + PC4 + (1|donor_id)
  
  # fit the random-effects variance model
  vp <- fitExtractVarPartModel(vobj, form_random, meta_sub)
  
  # Ensure UTF-8 is used (avoids ISOLatin1)
  options(encoding = "UTF-8")
  Sys.setlocale("LC_CTYPE", "en_US.UTF-8")
  
  # visualize variance fractions for top 10 genes
  pdf(file.path(outdir, paste0(ct, "_VarPart_top10.pdf")), width=8, height=6)
  plotVarPart(sortCols(vp[1:10, ]))
  dev.off()
  
  # violin plot: contribution of each variable across all genes
  pdf(file.path(outdir, paste0(ct, "_VarPart_violin.pdf")), width=8, height=6)
  plotVarPart(vp)
  dev.off()

  
  # If I want to remove the residuals from the plot
  # Suppose 'vp' is your variancePartition result (a data.frame)
  head(vp)
  
  # Drop the "Residuals" column
  vp_noresid <- vp[, colnames(vp) != "Residuals"]
  
  # Renormalize so that the fractions sum to 1 across remaining covariates
  vp_noresid <- vp_noresid / rowSums(vp_noresid)
  
  # Plot again without residuals
  pdf(file.path(outdir, paste0(ct, "_VarPart_violin_noResiduals.pdf")), width=8, height=6)
  variancePartition::plotVarPart(vp_noresid) +
    ggtitle("Variance partition without residuals (normalized)")
  dev.off()
  
  # Same but only for top 10 genes
  # visualize variance fractions for top 10 genes
  pdf(file.path(outdir, paste0(ct, "_VarPart_top10_noResiduals.pdf")), width=8, height=6)
  plotVarPart(sortCols(vp_noresid[1:10, ]))
  dev.off()

  
  # (B) Expression of a gene by covariate
  # png(file.path("plots", paste0(ct, "_geneExpressionByAge.png")), width = 800, height = 600)
  # plotExpression(vobj, "TP53", meta_sub, "age")
  # dev.off()
  
  
  #### Gene filtering summary (before vs after filtering)
  n_before <- nrow(filtered_counts_list_fixed[[ct]])
  n_after <- nrow(vobj$E)
  
  pdf(file.path(outdir, paste0("gene_filtering_", ct, ".pdf")), width = 6, height = 5)
  barplot(
    c(Before = n_before, After = n_after),
    col = c("grey70", "lightgreen"),
    main = paste("Gene filtering -", ct),
    ylab = "Number of genes"
  )
  dev.off()
  
  message(" Completed ", ct)

  
  ### Plot model residuals to see whetherassumptions hold
  # If residuals are symetric and centerd at 0, model fits reasonably well
  pdf(file.path(outdir, paste0(ct, "_ResidualDistrib.pdf")), width=8, height=6)
  res <- residuals(fit_full)
  hist(res[1, ], breaks = 50, main = "Residual distribution for Gene 1")
  dev.off()


  ### Plot correlation  matx between pairs of covariates
  # meta_sub: your metadata for a single cell type
  
  bio_covs <- c("ancestry", "age_at_recruitment_first", "sex", "BMI")
  
  pdf(file.path(outdir, paste0(ct, "_CorrMatrix.pdf")), width=8, height=6)
  # cc <- canCorPairs(meta_sub[, bio_covs])
  form_random <- as.formula(paste("~", paste(bio_covs, collapse = " + ")))
  cc <- variancePartition::canCorPairs(form_random, data = meta_sub)
  plotCorrMatrix(cc)
  dev.off()
  
  
### Volcano plot per covariate
ct <- "B_naive"

# Load your saved model fit (.rds file)
fit_full <- readRDS("/home/ivm/LMM_DGE_results/DGE_results_LRT_by_drop//B_naive_fit_full.rds")  

# Confirm available covariates
covariates <- colnames(fit_full$coefficients)
print(covariates)

# === Loop through each covariate and create a volcano plot ===
# === Volcano plots per covariate ===
for (cov in covariates) {
  if (cov %in% c("(Intercept)")) next  # skip intercept
  
  df <- data.frame(
    gene = rownames(fit_full$coefficients),
    logFC = fit_full$coefficients[, cov],
    P.Value = fit_full$p.value[, cov]
  ) %>%
    mutate(
      FDR = p.adjust(P.Value, method = "fdr"),
      negLogP = -log10(P.Value),
      sigFDR = FDR < 0.05
    )
  
  # Label top genes (you can adjust n)
  top_labels <- df %>%
    filter(sigFDR) %>%
    arrange(FDR) %>%
    slice_head(n = 15)  # label top 15 most significant genes
  
  p <- ggplot(df, aes(x = logFC, y = negLogP)) +
    geom_point(aes(color = sigFDR), alpha = 0.6, size = 1.2) +
    geom_text_repel(
      data = top_labels,
      aes(label = gene),
      size = 3,
      max.overlaps = 15,
      box.padding = 0.3,
      point.padding = 0.2,
      segment.color = "grey50"
    ) +
    scale_color_manual(values = c("grey70", "red")) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black") +
    labs(
      title = paste0("DGE by ", cov, ct),
      x = "log2 Fold Change",
      y = "-log10(p-value)",
      color = "FDR < 0.05"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top")
  
  # Save plot as PDF
  pdf(file.path(outdir, paste0("Volcano_", cov, "_", ct, ".pdf")),
      width = 7, height = 6)
  print(p)
  dev.off()
  
  # Print quick summary
  cat("\n", cov, ":",
      sum(df$sigFDR, na.rm = TRUE), "genes at FDR < 0.05\n")
}









# === Volcano plots per covariate ===
# Confirm available covariates
covariates <- colnames(fit_full$coefficients)
print(covariates)

# === Loop through each covariate and create a volcano plot ===
# === Volcano plots per covariate ===
for (cov in covariates) {
  if (cov %in% c("(Intercept)")) next  # skip intercept
  
  # Extract logFC and p-values from fit_full
  # Assuming fit_full is an MArrayLM object from dream/limma
  logFC <- fit_full$coefficients[, cov]   # replace with your covariate
  pvals <- fit_full$p.value[, cov]
  
  # Adjust p-values (FDR)
  fdr <- p.adjust(pvals, method = "fdr")
  
  # Create a data frame for plotting
  volcano_df <- data.frame(
    Gene = rownames(fit_full$coefficients),
    logFC = logFC,
    FDR = fdr
  )
  
  # Define thresholds
  logfc_cutoff <- 0.5
  fdr_cutoff <- 0.05
  
  # Add significance column
  volcano_df$Significance <- "NotSig"
  volcano_df$Significance[volcano_df$FDR <= fdr_cutoff & volcano_df$logFC >= logfc_cutoff] <- "Up"
  volcano_df$Significance[volcano_df$FDR <= fdr_cutoff & volcano_df$logFC <= -logfc_cutoff] <- "Down"
  
  # Volcano plot
  pdf(file.path(outdir, paste0("Volcano_color_", cov, ct, ".pdf")),
      width = 7, height = 6)
  ggplot(volcano_df, aes(x = logFC, y = -log10(FDR), color = Significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NotSig" = "grey")) +
    geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed", color = "black") +
    labs(title = paste0("DGE by", cov, ct), x = "Log2 Fold Change", y = "-log10(FDR)") +
    theme_minimal(base_size = 14)
  
  dev.off()
  
}
 





































    




















