#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# WDP BYM INLA Data Validation Utilities
# Comprehensive data validation and quality control functions
# Author: WDP Analysis Team
# Date: 2024

#' Validate input data completeness and structure
#' @param data_list List of loaded data frames
#' @param config Configuration list
#' @return List with validation results
validate_input_data <- function(data_list, config) {
  validation_results <- list(
    valid = TRUE,
    errors = character(0),
    warnings = character(0),
    info = character(0)
  )
  
  # Check required data frames
  required_datasets <- c("pca", "cdc", "pesticide", "adjacency", "mapping")
  missing_datasets <- setdiff(required_datasets, names(data_list))
  
  if (length(missing_datasets) > 0) {
    validation_results$valid <- FALSE
    validation_results$errors <- c(
      validation_results$errors,
      sprintf("Missing required datasets: %s", paste(missing_datasets, collapse = ", "))
    )
  }
  
  # Validate each dataset structure
  for (dataset_name in names(data_list)) {
    dataset <- data_list[[dataset_name]]
    
    # Check if dataset is not empty
    if (nrow(dataset) == 0) {
      validation_results$valid <- FALSE
      validation_results$errors <- c(
        validation_results$errors,
        sprintf("Dataset '%s' is empty", dataset_name)
      )
    }
    
    # Dataset-specific validations
    if (dataset_name == "cdc") {
      required_cols <- c("COUNTY_FIPS", "Year", "Deaths", "Population")
      missing_cols <- setdiff(required_cols, names(dataset))
      if (length(missing_cols) > 0) {
        validation_results$errors <- c(
          validation_results$errors,
          sprintf("CDC data missing columns: %s", paste(missing_cols, collapse = ", "))
        )
      }
    }
    
    if (dataset_name == "pesticide") {
      if (!"COUNTY_FIPS" %in% names(dataset)) {
        validation_results$errors <- c(
          validation_results$errors,
          "Pesticide data missing COUNTY_FIPS column"
        )
      }
      
      if (!("Year" %in% names(dataset) || "YEAR" %in% names(dataset))) {
        validation_results$errors <- c(
          validation_results$errors,
          "Pesticide data missing Year/YEAR column"
        )
      }
    }
  }
  
  return(validation_results)
}

#' Check data completeness for analysis requirements
#' @param model_data Prepared model data
#' @param min_counties Minimum number of counties required
#' @param min_records Minimum number of records required
#' @return Boolean indicating if data meets requirements
check_data_completeness <- function(model_data, min_counties = 50, min_records = 100) {
  if (is.null(model_data) || nrow(model_data) == 0) {
    return(FALSE)
  }
  
  n_counties <- length(unique(model_data$COUNTY_FIPS))
  n_records <- nrow(model_data)
  
  return(n_counties >= min_counties && n_records >= min_records)
}

#' Validate model requirements before fitting
#' @param model_data Model data frame
#' @param formula Model formula
#' @param config Configuration list
#' @return List with validation results
validate_model_requirements <- function(model_data, formula, config) {
  validation_results <- list(
    valid = TRUE,
    errors = character(0),
    warnings = character(0)
  )
  
  # Check response variable
  if (!"Deaths" %in% names(model_data)) {
    validation_results$valid <- FALSE
    validation_results$errors <- c(
      validation_results$errors,
      "Response variable 'Deaths' not found in model data"
    )
  }
  
  # Check exposure variable
  if (!"pesticide_lagged" %in% names(model_data)) {
    validation_results$valid <- FALSE
    validation_results$errors <- c(
      validation_results$errors,
      "Exposure variable 'pesticide_lagged' not found in model data"
    )
  }
  
  # Check for missing values in key variables
  key_vars <- c("Deaths", "Population", "pesticide_lagged")
  for (var in key_vars) {
    if (var %in% names(model_data)) {
      na_count <- sum(is.na(model_data[[var]]))
      if (na_count > 0) {
        validation_results$warnings <- c(
          validation_results$warnings,
          sprintf("Variable '%s' has %d missing values", var, na_count)
        )
      }
    }
  }
  
  return(validation_results)
}

#' Validate result row structure
#' @param result_row Result data frame row
#' @param config Configuration list
#' @return List with validation results
validate_result_row <- function(result_row, config) {
  validation_results <- list(
    valid = TRUE,
    errors = character(0),
    warnings = character(0)
  )
  
  # Check required columns
  required_cols <- c(
    "Timestamp", "Disease", "Exposure", "Category", "Measure", 
    "Estimate", "Lag", "Model", "Q1", "Q2", "Q3", "Q4", "Q5", 
    "P_Value", "N_Counties", "N_Records"
  )
  
  missing_cols <- setdiff(required_cols, names(result_row))
  if (length(missing_cols) > 0) {
    validation_results$valid <- FALSE
    validation_results$errors <- c(
      validation_results$errors,
      sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", "))
    )
  }
  
  # Validate specific values
  if ("P_Value" %in% names(result_row) && !is.na(result_row$P_Value)) {
    p_numeric <- suppressWarnings(as.numeric(gsub("\\*+", "", result_row$P_Value)))
    if (!is.na(p_numeric) && (p_numeric < 0 || p_numeric > 1)) {
      validation_results$warnings <- c(
        validation_results$warnings,
        sprintf("P-value outside valid range [0,1]: %.4f", p_numeric)
      )
    }
  }
  
  return(validation_results)
}

#' Print validation results
#' @param validation_results Results from validation functions
#' @param context Description of what was validated
print_validation_results <- function(validation_results, context = "Data") {
  cat(sprintf("📋 %s Validation Results:\n", context))
  
  if (validation_results$valid) {
    cat("  ✅ All validations passed\n")
  } else {
    cat("  ❌ Validation failed\n")
  }
  
  if (length(validation_results$errors) > 0) {
    cat("  🚨 Errors:\n")
    for (error in validation_results$errors) {
      cat(sprintf("    - %s\n", error))
    }
  }
  
  if (length(validation_results$warnings) > 0) {
    cat("  ⚠️  Warnings:\n")
    for (warning in validation_results$warnings) {
      cat(sprintf("    - %s\n", warning))
    }
  }
  
  if (length(validation_results$info) > 0) {
    cat("  ℹ️  Information:\n")
    for (info in validation_results$info) {
      cat(sprintf("    - %s\n", info))
    }
  }
  
  cat("\n")
}