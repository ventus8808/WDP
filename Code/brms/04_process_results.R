#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(purrr)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

read_config <- function() {
  this <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA)
  if (is.na(this)) {
    args <- commandArgs(trailingOnly = FALSE)
    this <- normalizePath(sub("^--file=", "", args[grep("--file=", args)][1]))
  }
  root <- dirname(dirname(dirname(this)))
  yaml::read_yaml(file.path(root, "config.yaml"))
}

main <- function() {
  cfg <- read_config()
  res <- cfg$brms_analysis$results
  out_dir <- res$csv_outputs

  files <- list.files(out_dir, pattern = "_fixef\\.csv$", full.names = TRUE)
  if (length(files) == 0) {
    message("No brms outputs to process.")
    return(invisible(NULL))
  }

  all_fe <- purrr::map_dfr(files, readr::read_csv, show_col_types = FALSE, .id = "src")
  all_fe <- all_fe %>%
    mutate(
      scenario = stringr::str_match(basename(.data$src), "_(.*)_fixef\\.csv$")[,2]
    ) %>%
    select(scenario, dplyr::everything(), -src)

  # Write combined fixed effects table
  comb_path <- file.path(res$reports, "brms_fixed_effects_combined.csv")
  dir.create(res$reports, showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(all_fe, comb_path)

  # Simple forest-like plot for main effects (excluding Intercept)
  p <- all_fe %>%
    filter(.data$term != "Intercept") %>%
    ggplot(aes(x = .data$Estimate, y = .data$term, xmin = .data$`Q2.5`, xmax = .data$`Q97.5`)) +
    geom_point() +
    geom_errorbarh(height = 0.2) +
  facet_wrap(~ scenario, scales = "free_y") +
    theme_minimal() +
    labs(title = "brms fixed effects (95% CrI)", x = "Effect", y = "Term")

  fig_path <- file.path(cfg$brms_analysis$results$figures, "brms_fixed_effects_forest.png")
  ggsave(fig_path, p, width = 10, height = 6, dpi = 150)
  message(sprintf("Wrote combined table: %s and figure: %s", comb_path, fig_path))
}

if (identical(environment(), globalenv())) {
  main()
}
