#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# WDP BYM INLA Data Processing Utilities
# Modular functions for data loading, processing, and preparation
# Author: WDP Analysis Team
# Date: 2024

#' Load and validate all required data files
#' @param config Configuration list from YAML
#' @param disease_code Disease code to analyze
#' @param measure_type Type of exposure measure ('Weight' or 'Density')
#' @return List of loaded data frames
load_all_data <- function(config, disease_code, measure_type = "Weight") {
  cat(sprintf("Loading all data files for disease: %s (Measure: %s)\n", disease_code, measure_type))

  tryCatch({
    # Construct file paths
    base_dir <- config$data_paths$base_dir

    # Choose pesticide data file based on measure_type
    pesticide_file <- if (measure_type == 'Density') {
      config$data_paths$pesticide_density_data
    } else {
      config$data_paths$pesticide_data
    }

    if (is.null(pesticide_file)) {
        stop(sprintf("Pesticide data path for measure_type '%s' not found in config.", measure_type))
    }

    paths <- list(
      pca = file.path(base_dir, config$data_paths$pca_covariates),
      cdc = file.path(base_dir, gsub("\\{disease_code\\}", disease_code, config$data_paths$cdc_data_template)),
      pesticide = file.path(base_dir, pesticide_file),
      adjacency = file.path(base_dir, config$data_paths$adjacency_data),
      mapping = file.path(base_dir, config$data_paths$pesticide_mapping)
    )

    # Validate file existence
    for (name in names(paths)) {
      if (!file.exists(paths[[name]])) {
        stop(sprintf("Required file not found: %s", paths[[name]]))
      }
    }

    # Load data with suppressed messages
    data_list <- suppressWarnings(suppressMessages({
      lapply(paths, function(p) {
        df <- read_csv(p, show_col_types = FALSE)
        cat(sprintf("  ✓ Loaded %s: %d records\n", basename(p), nrow(df)))
        return(df)
      })
    }))

    # Standardize year column name in pesticide data
    if ("YEAR" %in% names(data_list$pesticide)) {
      data_list$pesticide <- data_list$pesticide %>% rename(Year = YEAR)
    }

    # Validate data structure
    validate_data_structure(data_list, config)

    return(data_list)

  }, error = function(e) {
    cat(sprintf("✗ Error loading data: %s\n", e$message))
    return(NULL)
  })
}

#' Validate loaded data structure and completeness
#' @param data_list List of loaded data frames
#' @param config Configuration list
validate_data_structure <- function(data_list, config) {
  # Check PCA data has required columns
  pca_required <- c("COUNTY_FIPS", "Year", "SVI_PCA", "Climate_Factor_1", "Climate_Factor_2")
  missing_pca <- setdiff(pca_required, names(data_list$pca))
  if (length(missing_pca) > 0) {
    warning(sprintf("Missing PCA columns: %s", paste(missing_pca, collapse = ", ")))
  }

  # Check CDC data structure
  cdc_required <- c("COUNTY_FIPS", "Year", "Deaths", "Population")
  missing_cdc <- setdiff(cdc_required, names(data_list$cdc))
  if (length(missing_cdc) > 0) {
    stop(sprintf("Missing CDC columns: %s", paste(missing_cdc, collapse = ", ")))
  }

  # Check pesticide data structure for FIPS and a year column
  pesticide_has_fips <- "COUNTY_FIPS" %in% names(data_list$pesticide)
  pesticide_has_year <- "Year" %in% names(data_list$pesticide) || "YEAR" %in% names(data_list$pesticide)

  if (!pesticide_has_fips) {
    stop("Missing pesticide column: COUNTY_FIPS")
  }
  if (!pesticide_has_year) {
    stop("Missing pesticide column: Year or YEAR")
  }

  cat("✓ Data structure validation completed\n")
}

#' Manual rolling mean function (replaces zoo::rollmean)
#' @param x Numeric vector
#' @param k Window size
#' @param align Alignment ("right", "center", "left")
#' @return Numeric vector with rolling means
rolling_mean <- function(x, k, align = "right") {
  n <- length(x)
  result <- rep(NA, n)

  if (align == "right") {
    for (i in k:n) {
      result[i] <- mean(x[(i-k+1):i], na.rm = TRUE)
    }
  } else if (align == "center") {
    half_k <- floor(k/2)
    for (i in (half_k + 1):(n - half_k)) {
      result[i] <- mean(x[(i - half_k):(i + half_k)], na.rm = TRUE)
    }
  } else if (align == "left") {
    for (i in 1:(n-k+1)) {
      result[i] <- mean(x[i:(i+k-1)], na.rm = TRUE)
    }
  }

  return(result)
}

#' Calculate lagged exposure with moving averages
#' @param pesticide_data Pesticide data frame
#' @param pesticide_col_name Column name for specific pesticide
#' @param lag_years Number of lag years
#' @param config Configuration list
#' @return Data frame with lagged exposure
calculate_lagged_exposure <- function(pesticide_data, pesticide_col_name, lag_years, config) {
  cat(sprintf("  📊 Calculating %d-year lagged exposure for %s...\n", lag_years, pesticide_col_name))

  suppressWarnings({
    # Standardize column names
    year_col <- if("YEAR" %in% names(pesticide_data)) "YEAR" else "Year"
    pesticide_data_renamed <- pesticide_data %>%
      rename(Year = !!sym(year_col))

    # Calculate lagged exposure
    lagged_data <- pesticide_data_renamed %>%
      arrange(COUNTY_FIPS, Year) %>%
      group_by(COUNTY_FIPS) %>%
      mutate(
        pesticide_lagged = rolling_mean(.data[[pesticide_col_name]], lag_years, align = "right")
      ) %>%
      filter(Year >= (config$data_processing$year_range$start + lag_years - 1)) %>%
      ungroup()
  })

  cat(sprintf("  ✓ Lagged exposure calculated: %d records\n", nrow(lagged_data)))
  return(lagged_data)
}

#' Create exposure groups based on quantiles (P25, P50, P75)
#' @param exposure_data Numeric vector of exposure values
#' @param config Configuration list
#' @return Factor with levels: low, medium, high
create_exposure_groups <- function(exposure_data, config) {
  # Remove missing values for quantile calculation
  clean_exposure <- exposure_data[!is.na(exposure_data)]

  if (length(clean_exposure) == 0) {
    warning("No valid exposure data for grouping")
    return(factor(rep("medium", length(exposure_data)),
                  levels = c("medium", "low", "high")))
  }

  # Calculate quantiles
  q_low <- quantile(clean_exposure, config$exposure_grouping$low_quantile, na.rm = TRUE)
  q_high <- quantile(clean_exposure, config$exposure_grouping$high_quantile, na.rm = TRUE)

  # Create groups
  exposure_group <- case_when(
    is.na(exposure_data) ~ NA_character_,
    exposure_data <= q_low ~ "low",
    exposure_data >= q_high ~ "high",
    TRUE ~ "medium"
  )

  # Convert to factor with medium as reference
  factor(exposure_group, levels = c("medium", "low", "high"))
}

#' Apply winsorization to exposure data
#' @param x Numeric vector
#' @param config Configuration list
#' @return Winsorized numeric vector
#' @deprecated This function has been replaced by log transformation in prepare_model_data
# apply_winsorization <- function(x, config) {
#   if (!config$data_processing$winsorization$enabled) {
#     return(x)
#   }
#
#   # Remove NA values for quantile calculation
#   clean_x <- x[!is.na(x)]
#
#   if (length(clean_x) == 0) {
#     return(x)
#   }
#
#   q_low <- quantile(clean_x, config$data_processing$winsorization$lower_quantile, na.rm = TRUE)
#   q_high <- quantile(clean_x, config$data_processing$winsorization$upper_quantile, na.rm = TRUE)
#
#   # Apply winsorization
#   x[x < q_low] <- q_low
#   x[x > q_high] <- q_high
#
#   # Additional check for extreme outliers (more than 10 times the 99th percentile)
#   # This helps with very extreme values that might still exist after winsorization
#   p99 <- quantile(clean_x, 0.99, na.rm = TRUE)
#   extreme_threshold <- p99 * 10
#
#   if (is.finite(extreme_threshold) && extreme_threshold > 0) {
#     x[x > extreme_threshold] <- extreme_threshold
#   }
#
#   # Additional check for extremely high values (more than 100 times the 95th percentile)
#   # This is specifically for categories like category 19 that have severe outliers
#   p95 <- quantile(clean_x, 0.95, na.rm = TRUE)
#   extreme_threshold_95 <- p95 * 100
#
#   if (is.finite(extreme_threshold_95) && extreme_threshold_95 > 0) {
#     x[x > extreme_threshold_95] <- extreme_threshold_95
#   }
#
#   return(x)
# }

#' Prepare model data for INLA analysis
#' @param data_list List of loaded data
#' @param pesticide_col_name Pesticide column name
#' @param lag_years Lag period
#' @param config Configuration list
#' @return Prepared model data frame
prepare_model_data <- function(data_list, pesticide_col_name, lag_years, config) {
  cat("  🔧 Preparing model data...\n")

  suppressWarnings({
    # Calculate lagged exposure
    pesticide_lagged <- calculate_lagged_exposure(
      data_list$pesticide,
      pesticide_col_name,
      lag_years,
      config
    )

    # Prepare base model data
    model_data <- data_list$cdc %>%
      filter(
        Year >= (config$data_processing$year_range$start + lag_years - 1),
        Year <= config$data_processing$year_range$end
      ) %>%
      group_by(COUNTY_FIPS, Year) %>%
      summarise(
        Deaths = sum(Deaths, na.rm = TRUE),
        Population = sum(Population, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      left_join(data_list$pca, by = c("COUNTY_FIPS", "Year")) %>%
      left_join(pesticide_lagged, by = c("COUNTY_FIPS", "Year"))

    # Calculate expected deaths (indirect standardization)
    model_data <- model_data %>%
      group_by(Year) %>%
      mutate(
        national_rate = sum(Deaths, na.rm = TRUE) / sum(Population, na.rm = TRUE),
        expected_deaths = Population * national_rate,
        log_expected = log(expected_deaths + 1e-6)
      ) %>%
      ungroup()

    # Apply log transformation to exposure data
    # Step 1: Determine a small positive constant c for handling zero values
    non_zero_min <- min(model_data$pesticide_lagged[model_data$pesticide_lagged > 0], na.rm = TRUE)
    c_constant <- if (is.finite(non_zero_min)) non_zero_min / 2 else 0.001

    # Step 2: Create log-transformed column
    model_data$pesticide_log <- log(model_data$pesticide_lagged + c_constant)

    # Step 3: Standardize the log-transformed variable
    model_data$pesticide_log_std <- as.numeric(scale(model_data$pesticide_log))

    # Step 3b: Create binned index for non-linear modeling
    # Discretize the continuous exposure into bins for rw2 model
    n_bins <- if (!is.null(config$model_fitting$nonlinear$n_bins)) {
      config$model_fitting$nonlinear$n_bins
    } else {
      20  # Default to 20 bins
    }

    # Create equal-width bins on the log scale
    breaks <- seq(min(model_data$pesticide_log_std, na.rm = TRUE),
                  max(model_data$pesticide_log_std, na.rm = TRUE),
                  length.out = n_bins + 1)

    # Add small buffer to ensure all values are included
    breaks[1] <- breaks[1] - 0.001
    breaks[length(breaks)] <- breaks[length(breaks)] + 0.001

    # Create binned index
    model_data$pesticide_binned_idx <- as.numeric(cut(model_data$pesticide_log_std,
                                                       breaks = breaks,
                                                       include.lowest = TRUE))

    # Step 4: Robust outlier identification strategy
    # Calculate IQR-based outlier bounds on log-transformed data
    Q1 <- quantile(model_data$pesticide_log, 0.25, na.rm = TRUE)
    Q3 <- quantile(model_data$pesticide_log, 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1

    # Define outlier thresholds (3 times IQR)
    upper_bound <- Q3 + 3 * IQR
    lower_bound <- Q1 - 3 * IQR

    # Create outlier indicator
    model_data$is_outlier <- ifelse(
      model_data$pesticide_log > upper_bound | model_data$pesticide_log < lower_bound,
      TRUE, FALSE
    )

    # Report outlier statistics
    n_outliers <- sum(model_data$is_outlier, na.rm = TRUE)
    if (n_outliers > 0) {
      cat(sprintf("  ⚠️  Identified %d potential outliers (%.1f%%) based on IQR method\n",
                  n_outliers, 100 * n_outliers / nrow(model_data)))
    }

    # For compatibility with existing code, create pesticide_lagged_clean
    # This will be the non-outlier values from the original scale
    model_data$pesticide_lagged_clean <- ifelse(
      model_data$is_outlier,
      NA,
      model_data$pesticide_lagged
    )

    # Create exposure groups for backward compatibility and sensitivity analysis
    # Note: Main analysis will use continuous pesticide_log_std variable
    model_data$exposure_group <- create_exposure_groups(
      model_data$pesticide_lagged_clean,
      config
    )

    # Standardize covariates
    model_data <- model_data %>%
      mutate(
        # Standardize PCA variables (handle different naming conventions)
        SVI_std = if("SVI_PC1" %in% names(.)) as.numeric(scale(SVI_PC1)) else
                 if("SVI_PCA" %in% names(.)) as.numeric(scale(SVI_PCA)) else NA,
        Climate1_std = if("ENV_PC1" %in% names(.)) as.numeric(scale(ENV_PC1)) else
                      if("Climate_Factor_1" %in% names(.)) as.numeric(scale(Climate_Factor_1)) else NA,
        Climate2_std = if("ENV_PC2" %in% names(.)) as.numeric(scale(ENV_PC2)) else
                      if("Climate_Factor_2" %in% names(.)) as.numeric(scale(Climate_Factor_2)) else NA
      ) %>%
      filter(
        !is.na(exposure_group),
        !is.na(log_expected),
        is.finite(log_expected)
      )

    # Add county indices for spatial modeling
    counties <- sort(unique(model_data$COUNTY_FIPS))
    county_map <- data.frame(
      COUNTY_FIPS = counties,
      county_idx = 1:length(counties)
    )
    model_data <- model_data %>% left_join(county_map, by = "COUNTY_FIPS")

    # Quality control checks
    validate_model_data(model_data, config)

    cat(sprintf("  ✓ Model data prepared: %d records across %d counties\n",
                nrow(model_data), length(counties)))

    return(model_data)
  })
}

#' Validate prepared model data
#' @param model_data Prepared model data frame
#' @param config Configuration list
validate_model_data <- function(model_data, config) {
  # Check minimum sample sizes
  if (nrow(model_data) < config$data_processing$min_thresholds$records_per_analysis) {
    warning(sprintf("Low record count: %d (minimum: %d)",
                    nrow(model_data),
                    config$data_processing$min_thresholds$records_per_analysis))
  }

  n_counties <- length(unique(model_data$COUNTY_FIPS))
  if (n_counties < config$data_processing$min_thresholds$counties_per_analysis) {
    warning(sprintf("Low county count: %d (minimum: %d)",
                    n_counties,
                    config$data_processing$min_thresholds$counties_per_analysis))
  }

  # Check exposure group distribution
  group_counts <- table(model_data$exposure_group)
  cat(sprintf("  📊 Exposure groups: Low=%d, Medium=%d, High=%d\n",
              group_counts["low"], group_counts["medium"], group_counts["high"]))

  # Warn if any group is too small
  min_group_size <- min(group_counts)
  if (min_group_size < 50) {
    warning(sprintf("Small exposure group detected (n=%d). Results may be unstable.", min_group_size))
  }
}

#' Get pesticide categories from mapping data
#' @param mapping_data Pesticide mapping data frame
#' @param category_filter Filter for specific categories ("ALL", "TOP5", or specific names)
#' @return Filtered data frame with categories to analyze
get_pesticide_categories <- function(mapping_data, category_filter = "ALL") {
  categories <- mapping_data %>%
    select(cat_id = category1_id, category = category1_name) %>%
    distinct() %>%
    arrange(cat_id)

  # Convert cat_id to numeric if it's character
  categories$cat_id <- as.numeric(as.character(categories$cat_id))

  if (category_filter == "ALL") {
    return(categories)
  } else if (category_filter == "TOP5") {
    return(categories %>% head(5))
  } else if (category_filter == "TEST") {
    return(categories %>% head(1))
  } else {
    # Check if category_filter is numeric (cat_id) or text (category name)
    if (grepl("^\\d+$", category_filter)) {
      # Numeric category ID
      cat_id_filter <- as.numeric(category_filter)
      filtered <- categories %>% filter(cat_id == cat_id_filter)
    } else {
      # Text category name
      filtered <- categories %>% filter(category == category_filter)
    }

    if (nrow(filtered) == 0) {
      stop(sprintf("Category not found: %s", category_filter))
    }
    return(filtered)
  }
}

#' Get pesticide compounds from mapping data
#' @param mapping_data Pesticide mapping data frame
#' @param compound_filter Filter for specific compounds ("ALL" or specific compound IDs)
#' @return Filtered data frame with compounds to analyze
get_pesticide_compounds <- function(mapping_data, compound_filter = "ALL") {
  # For compound analysis, we need to work with compound-level data
  compounds <- mapping_data %>%
    select(compound_id, compound_name, cat_id = category1_id, category = category1_name) %>%
    arrange(compound_id)

  # Convert IDs to numeric if they're character
  compounds$compound_id <- as.numeric(as.character(compounds$compound_id))
  compounds$cat_id <- as.numeric(as.character(compounds$cat_id))

  if (compound_filter == "ALL") {
    return(compounds)
  } else if (compound_filter == "TEST") {
    return(compounds %>% head(1))
  } else {
    # Check if compound_filter is numeric (compound_id)
    if (grepl("^\\d+$", compound_filter)) {
      # Numeric compound ID
      compound_id_filter <- as.numeric(compound_filter)
      filtered <- compounds %>% filter(compound_id == compound_id_filter)
    } else {
      # Text compound name
      filtered <- compounds %>% filter(compound_name == compound_filter)
    }

    if (nrow(filtered) == 0) {
      stop(sprintf("Compound not found: %s", compound_filter))
    }
    return(filtered)
  }
}

#' Check if pesticide column exists in data
#' @param pesticide_data Pesticide data frame
#' @param pesticide_col Column name to check
#' @return Logical indicating existence
check_pesticide_column <- function(pesticide_data, pesticide_col) {
  exists <- pesticide_col %in% names(pesticide_data)

  if (!exists) {
    available_cols <- names(pesticide_data)[grepl("^cat\\d+_", names(pesticide_data))]
    cat(sprintf("  ✗ Column %s not found\n", pesticide_col))
    if (length(available_cols) > 0) {
      cat(sprintf("  Available columns: %s\n", paste(head(available_cols, 5), collapse = ", ")))
    }
  }

  return(exists)
}

#' Create summary statistics for exposure data
#' @param exposure_data Numeric exposure vector
#' @param exposure_name Name of the exposure
#' @return List of summary statistics
create_exposure_summary <- function(exposure_data, exposure_name) {
  clean_data <- exposure_data[!is.na(exposure_data)]

  if (length(clean_data) == 0) {
    return(list(
      name = exposure_name,
      n = 0,
      mean = NA,
      sd = NA,
      min = NA,
      max = NA,
      q25 = NA,
      q75 = NA
    ))
  }

  list(
    name = exposure_name,
    n = length(clean_data),
    mean = mean(clean_data),
    sd = sd(clean_data),
    min = min(clean_data),
    max = max(clean_data),
    q25 = quantile(clean_data, 0.25),
    q75 = quantile(clean_data, 0.75)
  )
}

#' Format numbers according to configuration
#' @param x Numeric value or vector
#' @param type Type of formatting ("rr", "p_value", "default")
#' @param config Configuration list
#' @return Formatted character vector
format_numbers <- function(x, type = "default", config = NULL) {
  if (is.null(config)) {
    # Default formatting
    if (type == "rr") {
      return(sprintf("%.4f", x))
    } else if (type == "p_value") {
      return(ifelse(x < 0.05, "<0.05", sprintf("%.2f", x)))
    } else {
      return(as.character(x))
    }
  }

  # Use config-based formatting
  if (type == "rr") {
    digits <- config$output$precision$rr
    return(sprintf(paste0("%.", digits, "f"), x))
  } else if (type == "p_value") {
    threshold <- config$output_format$p_value_threshold
    digits <- config$output$precision$p_value
    return(ifelse(x < threshold,
                  paste0("<", threshold),
                  sprintf(paste0("%.", digits, "f"), x)))
  } else {
    return(as.character(x))
  }
}

cat("✓ Data processing utilities loaded successfully\n")
