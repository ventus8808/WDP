#!/usr/bin/env Rscript
# Quick inspection of all ridgeline .rds files

files <- list.files("Result/Ridgeline", pattern = "Lag.*\\.rds$", full.names = TRUE)

for (f in files) {
  cat("\n========================================\n")
  cat("File:", basename(f), "\n")
  cat("========================================\n")
  
  data <- readRDS(f)
  
  cat("Cancer:     ", data$metadata$cancer_type, "\n")
  cat("EQI period: ", data$metadata$eqi_period, "\n")
  cat("AAMR period:", data$metadata$aamr_period, "\n")
  cat("Lag:        ", data$metadata$lag, "\n")
  cat("Model type: ", data$metadata$model_type, "\n")
  cat("N draws:    ", data$metadata$n_draws, "\n")
  cat("N obs:      ", data$metadata$n_obs, "\n")
  
  cat("\nConvergence:\n")
  cat("  Max R-hat:      ", sprintf("%.4f", data$metadata$convergence$max_rhat), "\n")
  cat("  Min ESS (bulk): ", sprintf("%.0f", data$metadata$convergence$min_ess_bulk), "\n")
  
  cat("\nSummary statistics:\n")
  print(data$summary[, c("quintile", "mean", "sd", "q025", "q975")])
  
  if (data$metadata$model_type == "multi_domain") {
    cat("\nDomains included:\n")
    print(table(data$draws_long$domain))
  }
}
