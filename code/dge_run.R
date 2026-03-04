#!/usr/bin/env Rscript

# ---- minimal argument parser (NO optparse dependency) ----
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("Missing value for ", flag)
  args[[i + 1]]
}

opt <- list(
  celltype         = get_arg("--celltype", "ALL"),
  main_r           = get_arg("--main_r", NULL),
  filtered_counts  = get_arg("--filtered_counts", NULL),
  pca_list         = get_arg("--pca_list", NULL),
  master_meta      = get_arg("--master_meta", NULL),
  phenos           = get_arg("--phenos", NULL)
)

required <- c("main_r","filtered_counts","pca_list","master_meta","phenos")
missing <- required[vapply(required, function(k) is.null(opt[[k]]) || opt[[k]] == "", logical(1))]
if (length(missing) > 0) {
  stop("Missing required args: ", paste(missing, collapse = ", "),
       "\nExample:\n  Rscript dge_run.R --celltype ALL --main_r code/scDEA_LM_CT3_UKB_batch.R ",
       "--filtered_counts filtered_counts_list.rds --pca_list pca_list.rds ",
       "--master_meta F3_UKB_adata_obs_with_metadata.csv --phenos cases_controls_all.tsv\n")
}

cat("Running DGE wrapper\n")
cat("  celltype:        ", opt$celltype, "\n")
cat("  main_r:          ", opt$main_r, "\n")
cat("  filtered_counts: ", opt$filtered_counts, "\n")
cat("  pca_list:        ", opt$pca_list, "\n")
cat("  master_meta:     ", opt$master_meta, "\n")
cat("  phenos:          ", opt$phenos, "\n")

# ---- Load inputs ----
filtered_counts_list <- readRDS(opt$filtered_counts)
pca_list             <- readRDS(opt$pca_list)
master_table_F3      <- read.csv(opt$master_meta, stringsAsFactors = FALSE)
phenotypes           <- read.delim(opt$phenos, stringsAsFactors = FALSE)

# ---- Put key objects into global env so your existing script can see them ----
assign("filtered_counts_list", filtered_counts_list, envir = .GlobalEnv)
assign("pca_list",            pca_list,            envir = .GlobalEnv)
assign("master_table_F3",     master_table_F3,     envir = .GlobalEnv)
assign("phenotypes",          phenotypes,          envir = .GlobalEnv)

# Optional: if your main script supports it
assign("CELLTYPE_TO_RUN", opt$celltype, envir = .GlobalEnv)

# ---- Source your pipeline script ----
source(opt$main_r)

cat("Completed DGE wrapper.\n")
