#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--main_r", type="character", default=NA,
              help="Path to main R script inside the repo, e.g. code/scDEA_LM_CT3_UKB.R")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.na(opt$main_r) || !nzchar(opt$main_r)) {
  stop("ERROR: --main_r is required")
}

cat("====================================================\n")
cat("Running main script:\n  ", opt$main_r, "\n")
cat("====================================================\n")

if (!file.exists(opt$main_r)) {
  stop("ERROR: main script not found: ", opt$main_r)
}

source(opt$main_r)

cat("\n====================================================\n")
cat("Main script completed successfully.\n")
cat("====================================================\n")
