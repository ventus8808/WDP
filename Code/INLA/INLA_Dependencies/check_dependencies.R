#!/usr/bin/env Rscript
# Quick dependency check for WDP INLA
# This script checks if all required R packages are installed

# Configure CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org/"))

required <- c("here", "dplyr", "readr", "yaml", "argparse", "progress", "INLA")
missing <- character(0)

cat("🔍 Checking required R packages...\n")

for (pkg in required) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    missing <- c(missing, pkg)
    cat(sprintf("  ❌ %s: NOT FOUND\n", pkg))
  } else {
    cat(sprintf("  ✅ %s: OK\n", pkg))
  }
}

if (length(missing) > 0) {
  cat("\n❌ Missing packages:", paste(missing, collapse = ", "), "\n")
  cat("📦 Run: Rscript INLA_Dependencies/install_packages.R\n")
  quit(status = 1)
} else {
  cat("\n✅ All dependencies satisfied\n")
  quit(status = 0)
}