#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# WDP BYM INLA Dashboard Printing Utilities
# Functions for creating a structured, dashboard-like console output
# Author: WDP Analysis Team
# Date: 2024

# Load required libraries
suppressMessages({
  library(stringr)
})

#' Prints a main section header
#' @param title The title of the section
print_section_header <- function(title) {
  cat(paste0("\n", str_pad("=", 70, "right", "="), "\n"))
  cat(paste0(" ", title, "\n"))
  cat(paste0(str_pad("=", 70, "right", "="), "\n\n"))
}

#' Prints the initial analysis setup and configuration block
#' @param args Parsed command-line arguments
#' @param config Loaded configuration list
#' @param compounds_to_analyze Data frame of compounds to be analyzed
#' @param total_combinations Total number of model combinations to run
print_analysis_setup <- function(args, config, compounds_to_analyze, total_combinations) {
  # Section 1: Setup & Configuration
  disease_name <- names(config$analysis$disease_codes)[sapply(config$analysis$disease_codes, function(x) x == args$disease_code)]
  title <- sprintf("WDP BYM INLA Production Analysis - %s (%s)", args$disease_code, disease_name)
  print_section_header(title)

  cat("[1] SETUP & CONFIGURATION\n")
  cat(sprintf("- Exposure Group: %s\n", args$pesticide_category))
  cat(sprintf("- Analysis Scope: %d combinations (%s measures × %s estimates × %s lags × %s models)\n",
              total_combinations,
              length(parse_list_argument(args$measure_type)),
              length(parse_list_argument(args$estimate_types)),
              length(parse_list_argument(args$lag_years)),
              length(parse_list_argument(args$model_types))))
  cat(sprintf("- Output Template: %s\n", file.path(config$output$base_dir, config$output$filename_template)))

  # Section 2: Data Sources
  cat("\n[2] DATA SOURCES\n")
  cat(sprintf("- Health Data: %s\n", file.path(config$data_paths$base_dir, gsub("\\{disease_code\\}", args$disease_code, config$data_paths$cdc_data_template))))
  cat(sprintf("- Exposure (Weight): %s\n", file.path(config$data_paths$base_dir, config$data_paths$pesticide_data)))
  cat(sprintf("- Exposure (Density): %s\n", file.path(config$data_paths$base_dir, config$data_paths$pesticide_density_data)))
  cat(sprintf("- Covariates: %s\n", file.path(config$data_paths$base_dir, config$data_paths$pca_covariates)))
  cat(sprintf("- Adjacency: %s\n", file.path(config$data_paths$base_dir, config$data_paths$adjacency_data)))
}

#' Prints a header for a specific compound being processed
#' @param index Current compound index
#' @param total Total number of compounds
#' @param compound_name The name of the compound
print_compound_header <- function(index, total, compound_name) {
  cat(paste0("\n", str_pad("-", 98, "right", "-"), "\n"))
  cat(sprintf("COMPOUND %d of %d: %s\n", index, total, compound_name))
}

#' Prints a summary for a completed compound analysis
#' @param output_path Path to the compound's result file
#' @param successful_runs Number of successful combinations
#' @param total_runs Total combinations for the compound
print_compound_summary <- function(output_path, successful_runs, total_runs) {
  success_rate <- (successful_runs / total_runs) * 100
  cat(sprintf("└─> Compound Complete. Success: %.1f%% (%d/%d). Results saved to: %s\n",
              success_rate, successful_runs, total_runs, output_path))
}

#' Prints the header for the real-time results table
print_results_table_header <- function() {
  cat("\n[3] REAL-TIME RESULTS & PROGRESS\n")
  header <- sprintf("| %-7s | %-8s | %-3s | %-5s | %-5s | %-30s | %-10s | %-8s | %-8s | %-7s |",
                    "Measure", "Estimate", "Lag", "Model", "N", "RR per SD [95% CI]", "P-Value", "DIC", "WAIC", "Status")
  cat(paste0(header, "\n"))
  cat(paste0(str_pad("-", nchar(header), "right", "-"), "\n"))
}

#' Formats and prints a single row of the real-time results table
#' @param result_row A single row data frame from create_result_row or create_failed_result_row
print_results_table_row <- function(result_row) {
  # Determine status by checking if RR_Per_SD is NA (a sign of failure)
  success <- !is.na(result_row$RR_Per_SD)
  status_icon <- if (success) "✅ OK" else "❌ ERROR"

  # Extract key information for display
  rr_str <- "NA"
  p_val_str <- "NA"

  if (success) {
    # Format the dose-response RR with confidence interval
    rr_str <- sprintf("%s [%s, %s]",
                     result_row$RR_Per_SD,
                     result_row$RR_Per_SD_Lower,
                     result_row$RR_Per_SD_Upper)
    p_val_str <- result_row$P_Value
  }

  # Format DIC - robustly
  dic_str <- "NA"
  if ("DIC" %in% names(result_row)) {
    dic_val <- suppressWarnings(as.numeric(result_row$DIC))
    # Ensure dic_val is a single, non-NA number before formatting
    if (length(dic_val) == 1 && !is.na(dic_val)) {
      dic_str <- sprintf("%.0f", dic_val)
    }
  }

  # Format WAIC - robustly
  waic_str <- "NA"
  if ("WAIC" %in% names(result_row)) {
    waic_val <- suppressWarnings(as.numeric(result_row$WAIC))
    # Ensure waic_val is a single, non-NA number before formatting
    if (length(waic_val) == 1 && !is.na(waic_val)) {
      waic_str <- sprintf("%.0f", waic_val)
    }
  }

  row_str <- sprintf("| %-7s | %-8s | %-3s | %-5s | %-5d | %-30s | %-10s | %-8s | %-8s | %-7s |",
                     result_row$Measure,
                     result_row$Estimate,
                     paste0(result_row$Lag, "y"),
                     result_row$Model,
                     result_row$N_Records,
                     rr_str,
                     p_val_str,
                     dic_str,
                     waic_str,
                     status_icon)
  cat(paste0(row_str, "\n"))
}

#' Prints the final summary of the entire analysis run.
#' This function displays total time, success rate, and details of any failed combinations.
#' @param failed_combinations A list of failed analysis result rows.
#' @param total_successful Total number of successful combinations.
#' @param total_combinations Total number of combinations attempted.
#' @param start_time The start time of the analysis (from Sys.time()).
print_final_summary <- function(failed_combinations, total_successful, total_combinations, start_time) {
  end_time <- Sys.time()
  total_time <- difftime(end_time, start_time, units = "mins")
  success_rate <- (total_successful / total_combinations) * 100

  cat("\n[4] FINAL SUMMARY\n")
  cat(sprintf("- Total Analysis Time: %.2f minutes\n", total_time))
  cat(sprintf("- Overall Success Rate: %.1f%% (%d/%d)\n",
              success_rate, total_successful, total_combinations))

  if (length(failed_combinations) > 0) {
    cat(sprintf("- Failed Combinations: %d\n", length(failed_combinations)))
    # Use a simple loop to print details of failed runs
    for (i in 1:length(failed_combinations)) {
      fail_row <- failed_combinations[[i]]
      # Basic printing, assuming fail_row is a data frame/list with expected names
      cat(sprintf("  %d. Measure=%s, Estimate=%s, Lag=%s, Model=%s, Error=%s\n",
                  i,
                  fail_row$Measure,
                  fail_row$Estimate,
                  fail_row$Lag,
                  fail_row$Model,
                  fail_row$Status_Message))
    }
  } else {
    cat("- All combinations completed successfully.\n")
  }
}

cat("✓ Dashboard printing utilities loaded successfully\n")
