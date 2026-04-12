#!/usr/bin/env Rscript
# Merge split sensitivity output files into final tables.
# Reads:  {prefix}_sensitivity_{EQI|Covar}_{GroupTag}.csv  (one per model group × cancer)
# Writes: {prefix}_sensitivity_EQI.csv
#         {prefix}_sensitivity_Covar.csv
# Row order: AAMR_Period → Lag → Sensitivity_Model (canonical order preserved)

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(readr)
})

option_list <- list(
  make_option(c("--output-dir"), type="character", default="Result/brms_Sensativity")
)
opt <- parse_args(OptionParser(option_list=option_list))

project_root <- normalizePath(".")
out_dir <- file.path(project_root, opt$`output-dir`)

PREFIXES      <- c("NDD", "Cancer")
MODEL_GROUPS  <- c("Baseline_Full", "Insurance", "Doctor", "Forest", "Monitoring")
TABLE_TYPES   <- c("EQI", "Covar")

# Canonical row order for Sensitivity_Model
MODEL_ORDER <- c("Baseline", "Full", "+Insurance", "+Doctor", "+Forest", "+Monitoring")

merge_ok <- TRUE

for (prefix in PREFIXES) {
  for (ttype in TABLE_TYPES) {
    split_files <- file.path(out_dir,
      paste0(prefix, "_sensitivity_", ttype, "_", MODEL_GROUPS, ".csv"))

    present <- split_files[file.exists(split_files)]
    missing <- split_files[!file.exists(split_files)]

    if (length(missing) > 0) {
      message("⚠️  Missing split files for ", prefix, " ", ttype, ":")
      for (f in missing) message("    ", basename(f))
      merge_ok <- FALSE
    }

    if (length(present) == 0) {
      message("⚠️  No split files found for ", prefix, " ", ttype, ", skipping merge.")
      next
    }

    combined <- rbindlist(lapply(present, fread), fill=TRUE)

    # Sort: AAMR_Period → Lag → Sensitivity_Model (canonical order)
    combined[, Sensitivity_Model := factor(Sensitivity_Model, levels=MODEL_ORDER)]
    combined <- combined[order(AAMR_Period, Lag, Sensitivity_Model)]
    combined[, Sensitivity_Model := as.character(Sensitivity_Model)]

    out_file <- file.path(out_dir, paste0(prefix, "_sensitivity_", ttype, ".csv"))
    write_csv(combined, out_file)
    message("✓ ", basename(out_file), "  (", nrow(combined), " rows from ",
            length(present), " split files)")
  }
}

if (merge_ok) {
  message("\nAll merges complete. Cleaning up split files...")
  for (prefix in PREFIXES) {
    for (ttype in TABLE_TYPES) {
      for (grp in MODEL_GROUPS) {
        f <- file.path(out_dir, paste0(prefix, "_sensitivity_", ttype, "_", grp, ".csv"))
        if (file.exists(f)) file.remove(f)
      }
    }
  }
  message("Split files removed.")
} else {
  message("\n⚠️  Some split files were missing — split files NOT removed. Re-run failed tasks then merge again.")
}
