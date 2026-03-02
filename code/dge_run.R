#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--celltype", type="character"),
  make_option("--main_r", type="character"),
  make_option("--filtered_counts", type="character"),
  make_option("--pca_list", type="character"),
  make_option("--master_meta", type="character"),
  make_option("--phenos", type="character")
)

opt <- parse_args(OptionParser(option_list = option_list))

cat("Running DGE for celltype:", opt$celltype, "\n")
cat("Main script:", opt$main_r, "\n")

# Load inputs
filtered_counts_list <- readRDS(opt$filtered_counts)
pca_list <- readRDS(opt$pca_list)

# Your master metadata + phenotypes (adapt to your actual reading functions)
master_table_F3 <- read.csv(opt$master_meta, stringsAsFactors = FALSE)
phenotypes <- read.delim(opt$phenos, stringsAsFactors = FALSE)

# Put key objects into global env so your existing script can see them
assign("filtered_counts_list", filtered_counts_list, envir = .GlobalEnv)
assign("pca_list", pca_list, envir = .GlobalEnv)
assign("master_table_F3", master_table_F3, envir = .GlobalEnv)
assign("phenotypes", phenotypes, envir = .GlobalEnv)

# Also pass the celltype you want to run
assign("CELLTYPE_TO_RUN", opt$celltype, envir = .GlobalEnv)

# Source your pipeline script
source(opt$main_r)

# IMPORTANT:
# Your main script must be written to:
# - either run only CELLTYPE_TO_RUN if it exists, OR
# - expose a function like run_dge_for_celltype(ct, ...)
# If your script currently loops over all ct, add a guard inside it.

cat("Completed DGE for celltype:", opt$celltype, "\n")
