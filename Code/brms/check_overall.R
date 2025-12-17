#!/usr/bin/env Rscript
files <- list.files("Result/Ridgeline", pattern = "Overall\\.rds$", full.names = TRUE)
for (f in files) {
  cat("\n", basename(f), "\n")
  data <- readRDS(f)
  cat("Lag:", data$metadata$lag, "| Rhat:", sprintf("%.4f", data$metadata$convergence$max_rhat), "\n")
  print(data$summary[, c("quintile", "mean", "q025", "q975")])
}
