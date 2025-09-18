#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# WDP BYM INLA Result Extraction Utilities
# Functions for extracting, formatting, and outputting analysis results
# Author: WDP Analysis Team
# Date: 2024

#' Extract coefficient results from INLA model
#' @param model INLA model object
#' @param coef_name Coefficient name to extract
#' @param config Configuration list
#' @return List with RR, confidence intervals, and p-value
extract_coefficient <- function(model, coef_name, config = NULL) {
  if (is.null(model) || is.null(model$summary.fixed)) {
    return(list(RR = NA, RR_Lower = NA, RR_Upper = NA, P_value = NA))
  }

  fixed_effects <- model$summary.fixed

  if (!coef_name %in% rownames(fixed_effects)) {
    return(list(RR = NA, RR_Lower = NA, RR_Upper = NA, P_value = NA))
  }

  coef_data <- fixed_effects[coef_name, ]

  # Calculate relative risk and confidence interval
  rr <- exp(coef_data$mean)
  rr_lower <- exp(coef_data$`0.025quant`)
  rr_upper <- exp(coef_data$`0.975quant`)

  # Calculate p-value (two-tailed test)
  p_value <- 2 * (1 - pnorm(abs(coef_data$mean / coef_data$sd)))

  return(list(
    RR = rr,
    RR_Lower = rr_lower,
    RR_Upper = rr_upper,
    P_value = p_value
  ))
}

#' Extract continuous dose-response effect
#' @param model INLA model object
#' @param config Configuration list
#' @return List with continuous exposure effect per SD increase
extract_exposure_effects <- function(model, config = NULL) {
  # Check if this is a non-linear model (has pesticide_binned_idx in random effects)
  if (!is.null(model) && !is.null(model$summary.random) &&
      "pesticide_binned_idx" %in% names(model$summary.random)) {

    # Non-linear model: extract from random effects
    rw2_summary <- model$summary.random$pesticide_binned_idx
    n_bins <- nrow(rw2_summary)

    # Calculate RR for each bin relative to first bin (reference)
    bin_rr <- exp(rw2_summary$mean - rw2_summary$mean[1])
    bin_rr_lower <- exp(rw2_summary$`0.025quant` - rw2_summary$mean[1])
    bin_rr_upper <- exp(rw2_summary$`0.975quant` - rw2_summary$mean[1])

    # Estimate P90 vs P10 by finding bins corresponding to these percentiles
    p10_bin <- ceiling(n_bins * 0.1)
    p90_bin <- floor(n_bins * 0.9)

    # Calculate RR for P90 vs P10
    rr_p90_vs_p10 <- exp(rw2_summary$mean[p90_bin] - rw2_summary$mean[p10_bin])
    rr_p90_vs_p10_lower <- exp(rw2_summary$`0.025quant`[p90_bin] - rw2_summary$`0.975quant`[p10_bin])
    rr_p90_vs_p10_upper <- exp(rw2_summary$`0.975quant`[p90_bin] - rw2_summary$`0.025quant`[p10_bin])

    # Calculate approximate "per SD" effect as average slope
    # This is less meaningful for non-linear models but included for consistency
    total_change <- rw2_summary$mean[n_bins] - rw2_summary$mean[1]
    rr_per_sd_approx <- exp(total_change / (n_bins / 3))  # Approximate 3 SDs across range

    # Calculate p-value by testing if any bin differs significantly from reference
    # Using conservative approach: minimum p-value across all bins
    z_scores <- abs(rw2_summary$mean) / rw2_summary$sd
    p_values <- 2 * (1 - pnorm(z_scores))
    p_value <- min(p_values[2:n_bins])  # Exclude reference bin

    return(list(
      is_nonlinear = TRUE,
      RR_per_SD = rr_per_sd_approx,
      RR_per_SD_Lower = NA,  # Not meaningful for non-linear
      RR_per_SD_Upper = NA,
      RR_P90_vs_P10_direct = rr_p90_vs_p10,
      RR_P90_vs_P10_Lower_direct = rr_p90_vs_p10_lower,
      RR_P90_vs_P10_Upper_direct = rr_p90_vs_p10_upper,
      P_value = p_value,
      bin_effects = data.frame(
        bin = 1:n_bins,
        RR = bin_rr,
        RR_Lower = bin_rr_lower,
        RR_Upper = bin_rr_upper
      )
    ))
  } else {
    # Linear model: extract from fixed effects
    continuous_effect <- extract_coefficient(model, "pesticide_log_std", config)

    return(list(
      is_nonlinear = FALSE,
      RR_per_SD = continuous_effect$RR,
      RR_per_SD_Lower = continuous_effect$RR_Lower,
      RR_per_SD_Upper = continuous_effect$RR_Upper,
      P_value = continuous_effect$P_value
    ))
  }
}

#' Extract spatial random effects summary
#' @param model INLA model object
#' @param config Configuration list
#' @return List with spatial effect statistics
extract_spatial_effects <- function(model, config = NULL) {
  if (is.null(model) || is.null(model$summary.random) ||
      !"county_idx" %in% names(model$summary.random)) {
    return(list(RR = NA, RR_Lower = NA, RR_Upper = NA))
  }

  spatial_summary <- model$summary.random$county_idx

  # Calculate overall spatial effect as geometric mean of county effects
  spatial_means <- spatial_summary$mean
  geometric_mean <- exp(mean(spatial_means, na.rm = TRUE))

  # Calculate confidence bounds using quantiles of county effects
  spatial_lower <- exp(quantile(spatial_summary$`0.025quant`, 0.5, na.rm = TRUE))
  spatial_upper <- exp(quantile(spatial_summary$`0.975quant`, 0.5, na.rm = TRUE))

  return(list(
    RR = geometric_mean,
    RR_Lower = spatial_lower,
    RR_Upper = spatial_upper
  ))
}

#' Extract temporal random effects summary
#' @param model INLA model object
#' @param config Configuration list
#' @return List with temporal effect statistics
extract_temporal_effects <- function(model, config = NULL) {
  if (is.null(model) || is.null(model$summary.random) ||
      !"Year" %in% names(model$summary.random)) {
    return(list(RR = NA, RR_Lower = NA, RR_Upper = NA))
  }

  temporal_summary <- model$summary.random$Year

  # Calculate temporal trend as geometric mean of year effects
  temporal_means <- temporal_summary$mean
  geometric_mean <- exp(mean(temporal_means, na.rm = TRUE))

  # Calculate confidence bounds
  temporal_lower <- exp(quantile(temporal_summary$`0.025quant`, 0.5, na.rm = TRUE))
  temporal_upper <- exp(quantile(temporal_summary$`0.975quant`, 0.5, na.rm = TRUE))

  return(list(
    RR = geometric_mean,
    RR_Lower = temporal_lower,
    RR_Upper = temporal_upper
  ))
}

#' Extract all covariate effects from model
#' @param model INLA model object
#' @param config Configuration list
#' @return List with all covariate effects
extract_covariate_effects <- function(model, config = NULL) {
  covariates <- list(
    SVI = extract_coefficient(model, "SVI_std", config),
    ENV_PC1 = extract_coefficient(model, "Climate1_std", config),
    ENV_PC2 = extract_coefficient(model, "Climate2_std", config),
    ENV_PC3 = NA  # Placeholder - Climate_Factor_3 not in current data
  )

  # Set ENV_PC3 to NA structure if not available
  if (is.na(covariates$ENV_PC3)) {
    covariates$ENV_PC3 <- list(RR = NA, RR_Lower = NA, RR_Upper = NA, P_value = NA)
  }

  return(covariates)
}

#' Format numeric values according to output specifications
#' @param x Numeric value or vector
#' @param type Type of formatting ("rr", "p_value", "default")
#' @param config Configuration list
#' @return Formatted character vector
format_output_numbers <- function(x, type = "default", config = NULL) {
  if (is.null(x) || all(is.na(x))) {
    return(NA_character_)
  }

  if (type == "rr") {
    # RR values to 4 decimal places
    return(sprintf("%.4f", x))
  } else if (type == "p_value") {
    # P-values: <0.05 if less than 0.05, otherwise 2 decimal places
    return(ifelse(x < 0.05, "<0.05", sprintf("%.2f", x)))
  } else {
    return(as.character(x))
  }
}

#' Create a single result row for output
#' @param analysis_info List with analysis metadata
#' @param model INLA model object
#' @param model_data Data frame used for analysis
#' @param config Configuration list
#' @return Data frame with single result row
create_result_row <- function(analysis_info, model, model_data, config) {

  # Extract all effects from model
  exposure_effects <- extract_exposure_effects(model, config)
  spatial_effects <- extract_spatial_effects(model, config)
  temporal_effects <- extract_temporal_effects(model, config)
  covariate_effects <- extract_covariate_effects(model, config)

  # Get model diagnostics
  model_diagnostics <- get_model_diagnostics(model, config)

  # Calculate sample sizes
  n_counties <- length(unique(model_data$COUNTY_FIPS))
  n_records <- nrow(model_data)

  # Calculate P90 vs P10 RR
  rr_p90_vs_p10 <- NA
  rr_p90_vs_p10_lower <- NA
  rr_p90_vs_p10_upper <- NA

  # Determine dose-response type
  dose_response_type <- "Linear"

  if (!is.null(exposure_effects$is_nonlinear) && exposure_effects$is_nonlinear) {
    # Non-linear model: use directly calculated P90 vs P10
    dose_response_type <- "Non-linear_RW2"
    rr_p90_vs_p10 <- exposure_effects$RR_P90_vs_P10_direct
    rr_p90_vs_p10_lower <- exposure_effects$RR_P90_vs_P10_Lower_direct
    rr_p90_vs_p10_upper <- exposure_effects$RR_P90_vs_P10_Upper_direct
  } else if (!is.na(exposure_effects$RR_per_SD)) {
    # Linear model: estimate from per-SD effect
    # P90 is approximately 1.28 SDs above P10
    rr_p90_vs_p10 <- exposure_effects$RR_per_SD^1.28
    rr_p90_vs_p10_lower <- exposure_effects$RR_per_SD_Lower^1.28
    rr_p90_vs_p10_upper <- exposure_effects$RR_per_SD_Upper^1.28
  }

  # Format p-value with stars
  p_value_formatted <- format_p_value_with_stars(exposure_effects$P_value)

  # Create result row with new format
  result_row <- data.frame(
    Timestamp = format(Sys.time(), config$output_format$timestamp_format),
    Disease = analysis_info$disease_code,
    Exposure = analysis_info$exposure_name,
    Category = analysis_info$category_name,
    Measure = analysis_info$measure_type,
    Estimate = analysis_info$estimate_type,
    Lag = analysis_info$lag_years,
    Model = analysis_info$model_type,

    # Dose-response model type
    Dose_Response_Type = dose_response_type,

    # Linear model results
    RR_Per_SD = format_output_numbers(exposure_effects$RR_per_SD, "rr", config),
    RR_Per_SD_Lower = format_output_numbers(exposure_effects$RR_per_SD_Lower, "rr", config),
    RR_Per_SD_Upper = format_output_numbers(exposure_effects$RR_per_SD_Upper, "rr", config),

    # P90 vs P10 comparison (estimated from linear model)
    RR_P90_vs_P10 = format_output_numbers(rr_p90_vs_p10, "rr", config),
    RR_P90_vs_P10_Lower = format_output_numbers(rr_p90_vs_p10_lower, "rr", config),
    RR_P90_vs_P10_Upper = format_output_numbers(rr_p90_vs_p10_upper, "rr", config),

    # P-value with stars
    P_Value = p_value_formatted,

    # Model fit statistics
    DIC = model_diagnostics$dic,
    WAIC = model_diagnostics$waic,

    # Sample size information
    N_Counties = n_counties,
    N_Records = n_records,

    # Status
    Status_Message = "SUCCESS",

    stringsAsFactors = FALSE
  )

  return(result_row)
}

#' Create a failed analysis result row
#' @param analysis_info List with analysis metadata
#' @param error_message Error message string
#' @param config Configuration list
#' @return Data frame with failed result row
create_failed_result_row <- function(analysis_info, error_message = "Analysis failed", config) {
  result_row <- data.frame(
    Timestamp = format(Sys.time(), config$output_format$timestamp_format),
    Disease = analysis_info$disease_code,
    Exposure = analysis_info$exposure_name,
    Category = analysis_info$category_name,
    Measure = analysis_info$measure_type,
    Estimate = analysis_info$estimate_type,
    Lag = analysis_info$lag_years,
    Model = analysis_info$model_type,

    # Dose-response model type
    Dose_Response_Type = "Linear",  # Failed models default to Linear

    # Linear model results as NA for failed runs
    RR_Per_SD = NA_character_,
    RR_Per_SD_Lower = NA_character_,
    RR_Per_SD_Upper = NA_character_,

    # P90 vs P10 comparison as NA
    RR_P90_vs_P10 = NA_character_,
    RR_P90_vs_P10_Lower = NA_character_,
    RR_P90_vs_P10_Upper = NA_character_,

    # P-value as NA for failed runs
    P_Value = NA_character_,

    # Model fit statistics as NA
    DIC = NA,
    WAIC = NA,

    # Zero sample sizes
    N_Counties = 0,
    N_Records = 0,

    # Error status
    Status_Message = as.character(error_message),

    stringsAsFactors = FALSE
  )

  return(result_row)
}

#' Initialize output file with proper headers
#' @param output_path Full path to output file
#' @param config Configuration list
initialize_output_file <- function(output_path, config) {
  output_dir <- dirname(output_path)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  if (!file.exists(output_path) || !config$output$append_mode) {
    # Create header row with new dose-response format
    header_cols <- c(
      "Timestamp", "Disease", "Exposure", "Category", "Measure", "Estimate", "Lag", "Model",
      "Dose_Response_Type", "RR_Per_SD", "RR_Per_SD_Lower", "RR_Per_SD_Upper",
      "RR_P90_vs_P10", "RR_P90_vs_P10_Lower", "RR_P90_vs_P10_Upper",
      "P_Value", "DIC", "WAIC", "N_Counties", "N_Records", "Status_Message"
    )

    # Write header to file
    writeLines(paste(header_cols, collapse = ","), output_path)
  }
}

#' Write result row to output file
#' @param result_row Data frame with result row
#' @param output_path Path to output file
#' @param config Configuration list
write_result_row <- function(result_row, output_path, config) {
  suppressMessages({
    write_csv(result_row, output_path, append = TRUE)
  })
}

#' Print a summary of the result row to console
#' @param result_row Data frame row with results
print_result_summary <- function(result_row) {
  # Check if this is a successful run by checking if RR_Per_SD is not NA
  if (!is.na(result_row$RR_Per_SD)) {
    # Format the dose-response results for display
    dose_response_str <- sprintf("RR per SD: %s [%s, %s]",
                                result_row$RR_Per_SD,
                                result_row$RR_Per_SD_Lower,
                                result_row$RR_Per_SD_Upper)

    cat(sprintf("✓ Analysis successful: %s\n", dose_response_str))
    cat(sprintf("  P-Value: %s, DIC: %.2f, WAIC: %.2f\n",
                result_row$P_Value,
                ifelse(is.na(result_row$DIC), 0, result_row$DIC),
                ifelse(is.na(result_row$WAIC), 0, result_row$WAIC)))
  } else {
    cat("❌ Analysis failed\n")
    if (!is.na(result_row$Status_Message)) {
      cat(sprintf("  Error: %s\n", result_row$Status_Message))
    }
  }
}

#' Validate result row against configuration requirements
#' @param result_row Data frame with result row
#' @param config Configuration list
#' @return List with validation results
validate_result_row <- function(result_row, config) {
  validation_results <- list(valid = TRUE, messages = character(0))

  # Check required columns
  required_cols <- config$result_validation$required_columns
  missing_cols <- setdiff(required_cols, names(result_row))

  if (length(missing_cols) > 0) {
    validation_results$valid <- FALSE
    validation_results$messages <- c(
      validation_results$messages,
      sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", "))
    )
  }

  # Check for reasonable RR values if specified
  if (!is.null(config$quality_control$validation$rr_bounds)) {
    # Check dose-response RR values
    dose_response_cols <- c("RR_Per_SD", "RR_P90_vs_P10")
    for (col in dose_response_cols) {
      if (col %in% names(result_row) && !is.na(result_row[[col]])) {
        rr_val <- suppressWarnings(as.numeric(result_row[[col]]))
        rr_bounds <- config$quality_control$validation$rr_bounds

        if (!is.na(rr_val) && (rr_val < rr_bounds[1] || rr_val > rr_bounds[2])) {
          validation_results$messages <- c(
            validation_results$messages,
            sprintf("%s value outside reasonable bounds: %.4f", col, rr_val)
          )
        }
      }
    }
  }

  # Check P-value bounds
  if (!is.na(result_row$P_Value)) {
    # Extract numeric part of p-value (remove stars)
    p_numeric <- suppressWarnings(as.numeric(gsub("\\*+", "", result_row$P_Value)))
    p_bounds <- config$quality_control$validation$p_value_bounds

    if (!is.na(p_numeric) && (p_numeric < p_bounds[1] || p_numeric > p_bounds[2])) {
      validation_results$messages <- c(
        validation_results$messages,
        sprintf("P-value outside valid bounds: %.4f", p_numeric)
      )
    }
  }

  return(validation_results)
}

#' Create summary statistics across multiple results
#' @param results_file Path to results CSV file
#' @param config Configuration list
#' @return Data frame with summary statistics
create_results_summary <- function(results_file, config = NULL) {
  if (!file.exists(results_file)) {
    warning(sprintf("Results file not found: %s", results_file))
    return(NULL)
  }

  results <- read_csv(results_file, show_col_types = FALSE)

  if (nrow(results) == 0) {
    warning("No results found in file")
    return(NULL)
  }

  # Convert RR and P columns to numeric for analysis
  results$RR_numeric <- suppressWarnings(as.numeric(results$RR))
  results$P_numeric <- ifelse(results$P == "<0.05", 0.049,
                              suppressWarnings(as.numeric(results$P)))

  # Create summary by category and model
  summary_stats <- results %>%
    group_by(Category, Model) %>%
    summarise(
      n_analyses = n(),
      n_significant = sum(P_numeric < 0.05, na.rm = TRUE),
      mean_rr = mean(RR_numeric, na.rm = TRUE),
      median_rr = median(RR_numeric, na.rm = TRUE),
      min_p = min(P_numeric, na.rm = TRUE),
      max_rr = max(RR_numeric, na.rm = TRUE),
      min_rr = min(RR_numeric, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(min_p)

  return(summary_stats)
}

#' Format p-value with significance stars
#' @param p_value Numeric p-value
#' @return Formatted string with stars
format_p_value_with_stars <- function(p_value) {
  if (is.na(p_value)) {
    return(NA_character_)
  }

  if (p_value < 0.001) {
    return(sprintf("%.3f***", p_value))
  } else if (p_value < 0.01) {
    return(sprintf("%.3f**", p_value))
  } else if (p_value < 0.05) {
    return(sprintf("%.3f*", p_value))
  } else {
    return(sprintf("%.3f", p_value))
  }
}

cat("✓ Result extraction utilities loaded successfully\n")
