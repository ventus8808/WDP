#!/usr/bin/env Rscript
#' 
#' BRMS Model Runner - Refactored for Intelligent Self-Configuration
#' 
#' 重构后的设计原则：
#' 1. 接收单一参数：scenario_name (通过 --scenario 标志)
#' 2. 自我配置：从 config.yaml 读取场景的完整定义
#' 3. 数据准备：内置所有必要的筛选和变量选择逻辑
#' 4. 智能输出：基于场景配置生成适当的输出文件名
#' 5. 移除LMM兼容性逻辑：由04_process_results.R负责格式转换
#' 
#' 用法: Rscript 02_run_brms_model.R --scenario LungCancer_TotalEQI_Lag5_AllRUCC
#'

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(brms)
  library(posterior)
  library(tibble)
  library(stringr)
})

# Load utilities - handle script execution context
this_file <- tryCatch({
  normalizePath(sys.frame(1)$ofile)
}, error = function(e) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grep("^--file=", args)]
  if (length(file_arg) > 0) {
    normalizePath(sub("^--file=", "", file_arg[1]))
  } else {
    stop("Cannot determine script location")
  }
})
source(file.path(dirname(this_file), "utils.R"))

# Utility operator
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Parse command line arguments for scenario name
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  # Look for --scenario flag
  scenario_idx <- which(args == "--scenario")
  if (length(scenario_idx) == 0) {
    stop("Usage: Rscript 02_run_brms_model.R --scenario <scenario_name>")
  }
  
  if (scenario_idx == length(args)) {
    stop("--scenario flag requires a scenario name")
  }
  
  scenario_name <- args[scenario_idx + 1]
  message(sprintf("[brms] Target scenario: %s", scenario_name))
  
  return(list(scenario_name = scenario_name))
}

#' Find scenario configuration by name
find_scenario <- function(config, scenario_name) {
  scenarios <- config$brms_analysis$scenarios
  
  for (scenario in scenarios) {
    if (identical(scenario$name, scenario_name)) {
      return(scenario)
    }
  }
  
  stop(sprintf("Scenario '%s' not found in config. Available scenarios: %s", 
               scenario_name, 
               paste(sapply(scenarios, function(s) s$name), collapse = ", ")))
}

#' Resolve EQI column name from config.yaml mapping
resolve_eqi_column <- function(domain, rucc_filter, eqi_column_map) {
  # Normalize domain key
  domain_key <- if (is.null(domain) || domain == '' || is.na(domain)) 'total' else as.character(domain)
  
  # Check if RUCC filter exists
  has_rucc_filter <- !is.null(rucc_filter) && length(rucc_filter) > 0 && !all(is.na(rucc_filter))
  
  # Choose the appropriate mapping
  map_key <- if (has_rucc_filter) "with_rucc" else "no_rucc"
  domain_map <- eqi_column_map[[map_key]]
  
  if (is.null(domain_map)) {
    stop(sprintf("EQI column mapping for '%s' not found in config", map_key))
  }
  
  if (!domain_key %in% names(domain_map)) {
    available_domains <- paste(names(domain_map), collapse = ", ")
    stop(sprintf("Unknown domain '%s'. Available domains: %s", domain_key, available_domains))
  }
  
  column_name <- domain_map[[domain_key]]
  
  message(sprintf("[brms] Using EQI column: %s (domain=%s, rucc_filter=%s)", 
                  column_name, domain_key, ifelse(is.null(rucc_filter), "none", toString(rucc_filter))))
  
  return(column_name)
}

#' Build priors from config.yaml specification
build_priors <- function(prior_key, prior_configs, data, has_re = FALSE) {
  if (!prior_key %in% names(prior_configs)) {
    available_keys <- paste(names(prior_configs), collapse = ", ")
    stop(sprintf("Prior key '%s' not found. Available: %s", prior_key, available_keys))
  }
  
  prior_config <- prior_configs[[prior_key]]
  
  # Extract values
  intercept_mean <- prior_config$intercept$fallback_mean %||% 100.0
  intercept_sd <- prior_config$intercept$fallback_sd %||% 50.0
  
  # Auto-center intercept if requested
  if (prior_config$intercept$strategy == "auto_center") {
    center_val <- tryCatch({
      median(c(data$AAMR_lo2, data$AAMR_hi2), na.rm = TRUE)
    }, error = function(e) {
      intercept_mean
    })
    if (is.finite(center_val)) intercept_mean <- center_val
  }
  
  b_sd <- prior_config$fixed_effects$class_b_sd %||% 5.0
  sigma_df <- prior_config$sigma$df %||% 3
  sigma_scale <- prior_config$sigma$scale %||% 30.0
  
  # Build prior list using set_prior for better compatibility
  priors <- c(
    brms::set_prior(sprintf("normal(%g, %g)", intercept_mean, intercept_sd), class = "Intercept"),
    brms::set_prior(sprintf("normal(0, %g)", b_sd), class = "b"),
    brms::set_prior(sprintf("student_t(%g, 0, %g)", sigma_df, sigma_scale), class = "sigma")
  )
  
  # Add random effects priors if needed
  if (has_re) {
    re_config <- prior_config$random_effects$state_sd
    re_df <- re_config$df %||% 3
    re_scale <- re_config$scale %||% 20.0
    priors <- c(priors, 
      brms::set_prior(sprintf("student_t(%g, 0, %g)", re_df, re_scale), class = "sd")
    )
  }
  
  return(priors)
}

#' Build model formula from config.yaml definition
build_formula <- function(formula_key, formula_configs, eqi_col_active, field_map) {
  if (!formula_key %in% names(formula_configs)) {
    available_keys <- paste(names(formula_configs), collapse = ", ")
    stop(sprintf("Formula key '%s' not found. Available: %s", formula_key, available_keys))
  }
  
  formula_template <- formula_configs[[formula_key]]
  
  # Replace placeholder variables with actual column names
  state_col <- field_map$state %||% "State"
  eqi_total_col <- field_map$eqi_total %||% "EQI"
  
  # Dynamic substitutions
  formula_rhs <- formula_template
  formula_rhs <- gsub("EQI_quintile", eqi_col_active, formula_rhs, fixed = TRUE)
  formula_rhs <- gsub("State", state_col, formula_rhs, fixed = TRUE)
  formula_rhs <- gsub("\\bEQI\\b", eqi_total_col, formula_rhs)  # Word boundary to avoid partial matches
  
  message(sprintf("[brms] Built formula RHS: %s", formula_rhs))
  
  return(formula_rhs)
}



main <- function() {
  message("[brms] Starting intelligent brms model runner (refactored)")
  
  # Parse command line arguments
  args_list <- parse_args()
  scenario_name <- args_list$scenario_name
  
  # Load configuration using centralized function
  cfg <- load_project_config()
  brms_config <- cfg$brms_analysis
  stopifnot(!is.null(brms_config))
  
  # Find scenario configuration
  scenario <- find_scenario(cfg, scenario_name)
  message(sprintf("[brms] Found scenario: %s", scenario$name))
  message(sprintf("[brms]   Cancer: %s, Domain: %s, Formula: %s, Family: %s", 
                  scenario$cancer_type, scenario$domain, scenario$formula_key, scenario$family))
  
  # Ensure output directories exist
  results_config <- brms_config$results
  for (dir_path in unlist(results_config)) {
    dir.create(dir_path, showWarnings = FALSE, recursive = TRUE)
  }
  
  # Load data
  data_file <- brms_config$data_file
  message(sprintf("[brms] Loading data from: %s", data_file))
  df <- readr::read_csv(data_file, show_col_types = FALSE)
  message(sprintf("[brms] Loaded %d rows, %d columns", nrow(df), ncol(df)))
  
  # Get field mappings
  field_map <- brms_config$field_map
  state_col <- field_map$state %||% 'State'
  cancer_col <- field_map$cancer_type %||% 'Cancer_Type'
  lo_col <- field_map$aamr_lower %||% 'AAMR_lower'
  hi_col <- field_map$aamr_upper %||% 'AAMR_upper'
  rucc_col <- field_map$rucc %||% 'RUCC'
  
  # Resolve active EQI column based on scenario configuration
  eqi_active_col <- resolve_eqi_column(scenario$domain, scenario$rucc_filter, brms_config$eqi_column_map)

  
  # Data type conversions and preprocessing
  message("[brms] Preprocessing data...")
  
  # Convert EQI column to factor with Q1-Q5 labels if numeric
  if (eqi_active_col %in% names(df)) {
    if (is.numeric(df[[eqi_active_col]]) || is.integer(df[[eqi_active_col]])) {
      df[[eqi_active_col]] <- factor(paste0('Q', df[[eqi_active_col]]), levels = paste0('Q', 1:5))
      message(sprintf("[brms] Converted %s to Q1-Q5 factor", eqi_active_col))
    }
  } else {
    stop(sprintf("Required EQI column '%s' not found in data", eqi_active_col))
  }
  
  # Convert categorical columns to factors
  df <- df %>% mutate(across(all_of(intersect(c(state_col, cancer_col), names(df))), as.factor))
  
  # Create standardized smoking rate if not exists
  smoking_col <- field_map$smoking_rate %||% 'Smoking_Rate'
  if (!'Smoking_Rate_std' %in% names(df)) {
    df <- df %>% mutate(Smoking_Rate_std = as.numeric(scale(.data[[smoking_col]])))
  }
  
  # Create censoring indicator if not exists
  if (!'cens_indicator' %in% names(df)) {
    df <- df %>% mutate(cens_indicator = ifelse(.data[[lo_col]] == .data[[hi_col]], 'none', 'interval'))
  }
  
  # === DATA FILTERING BASED ON SCENARIO CONFIGURATION ===
  message("[brms] Applying scenario-based filters...")
  
  # Start with full dataset
  df_filtered <- df
  
  # Filter by cancer type
  df_filtered <- df_filtered %>% filter(.data[[cancer_col]] == scenario$cancer_type)
  message(sprintf("[brms] After cancer filter (%s): %d rows", scenario$cancer_type, nrow(df_filtered)))
  
  # Filter by EQI period
  if (!is.null(scenario$eqi_period) && "EQI_Period" %in% names(df_filtered)) {
    target_eqi <- scenario$eqi_period
    df_filtered <- df_filtered %>% filter(as.character(.data$EQI_Period) == as.character(target_eqi))
    message(sprintf("[brms] After EQI period filter (%s): %d rows", target_eqi, nrow(df_filtered)))
  }
  
  # Filter by time period
  if (!is.null(scenario$time_period) && "Time_Period" %in% names(df_filtered)) {
    target_time <- scenario$time_period
    df_filtered <- df_filtered %>% filter(as.character(.data$Time_Period) == as.character(target_time))
    message(sprintf("[brms] After time period filter (%s): %d rows", target_time, nrow(df_filtered)))
  }
  
  # Filter by lag years
  if (!is.null(scenario$lag_years) && "Lag_Years" %in% names(df_filtered)) {
    target_lag <- scenario$lag_years
    df_filtered <- df_filtered %>% filter(as.numeric(.data$Lag_Years) == as.numeric(target_lag))
    message(sprintf("[brms] After lag years filter (%s): %d rows", target_lag, nrow(df_filtered)))
  }
  
  # Filter by RUCC (if specified)
  if (!is.null(scenario$rucc_filter) && rucc_col %in% names(df_filtered)) {
    rucc_values <- scenario$rucc_filter
    if (length(rucc_values) == 1) {
      df_filtered <- df_filtered %>% filter(as.numeric(.data[[rucc_col]]) == as.numeric(rucc_values))
      message(sprintf("[brms] After RUCC filter (%s): %d rows", rucc_values, nrow(df_filtered)))
    } else {
      df_filtered <- df_filtered %>% filter(.data[[rucc_col]] %in% rucc_values) 
      message(sprintf("[brms] After RUCC filter (%s): %d rows", paste(rucc_values, collapse=","), nrow(df_filtered)))
    }
  }
  
  # Final assignment
  df_sub <- df_filtered
  
  # Check if we have any data left after filtering
  if (nrow(df_sub) == 0) {
    message("[brms] ERROR: No data remains after filtering!")
    # Show available combinations for debugging
    available_combos <- df %>% 
      dplyr::count(.data[[cancer_col]], .data$EQI_Period, .data$Time_Period, .data$Lag_Years, name = "n") %>%
      dplyr::arrange(.data[[cancer_col]], .data$EQI_Period, .data$Time_Period, .data$Lag_Years)
    message("[brms] Available data combinations:")
    print(available_combos)
    stop("Filtered dataset is empty; check scenario configuration")
  }
  
  message(sprintf("[brms] Final dataset: %d rows after all filters", nrow(df_sub)))
  
  # === DATA VALIDATION AND CLEANING ===
  message("[brms] Validating and cleaning data...")
  
  # Ensure no NA in required fields
  req_cols <- c(lo_col, hi_col, smoking_col, 'Smoking_Rate_std', eqi_active_col, state_col)
  keep_cols <- intersect(req_cols, names(df_sub))
  
  rows_before <- nrow(df_sub)
  df_sub <- df_sub %>% filter(if_all(all_of(keep_cols), ~ !is.na(.x)))
  rows_after <- nrow(df_sub)
  
  if (rows_before != rows_after) {
    message(sprintf("[brms] Removed %d rows with missing values", rows_before - rows_after))
  }
  
  if (nrow(df_sub) == 0) {
    stop("No valid data remains after removing missing values")
  }
  
  # Fix any inverted intervals (lower > upper)
  lo_vec <- df_sub[[lo_col]]
  hi_vec <- df_sub[[hi_col]]
  inverted_count <- sum(lo_vec > hi_vec, na.rm = TRUE)
  
  if (inverted_count > 0) {
    message(sprintf("[brms] Fixed %d inverted intervals", inverted_count))
  }
  
  df_sub$AAMR_lo2 <- pmin(lo_vec, hi_vec)
  df_sub$AAMR_hi2 <- pmax(lo_vec, hi_vec)
  
  # Recompute censoring indicator based on cleaned bounds
  df_sub$cens2 <- factor(
    ifelse(df_sub$AAMR_lo2 == df_sub$AAMR_hi2, 'none', 'interval'),
    levels = c('none','interval')
  )
  
  # Drop unused EQI factor levels to avoid model singularities
  if (is.factor(df_sub[[eqi_active_col]])) {
    original_levels <- nlevels(df_sub[[eqi_active_col]])
    df_sub[[eqi_active_col]] <- droplevels(df_sub[[eqi_active_col]])
    new_levels <- nlevels(df_sub[[eqi_active_col]])
    
    if (original_levels != new_levels) {
      message(sprintf("[brms] Dropped %d unused EQI factor levels", original_levels - new_levels))
    }
  }
  
  message(sprintf("[brms] Data validation complete. Final dataset: %d rows", nrow(df_sub)))

  
  # === MODEL SPECIFICATION ===
  message("[brms] Building model specification...")
  
  # Build formula from config specification
  formula_key <- scenario$formula_key
  formula_rhs <- build_formula(formula_key, brms_config$formulas, eqi_active_col, field_map)
  
  message(sprintf("[brms] Using formula: %s", formula_key))
  
  # Build response specification
  use_midpoint <- isTRUE(scenario$use_midpoint)
  
  if (use_midpoint) {
    # Use midpoint of interval for simpler modeling
    df_sub$response_mid <- (df_sub$AAMR_lo2 + df_sub$AAMR_hi2) / 2
    lhs <- 'response_mid'
    message("[brms] Using interval midpoint as response")
  } else {
    # Use full interval censoring
    lhs <- sprintf("%s | cens(cens2, %s)", 'AAMR_lo2', 'AAMR_hi2')
    message("[brms] Using interval-censored response")
  }
  
  # Combine into full formula
  full_formula <- stats::as.formula(paste(lhs, "~", formula_rhs))
  
  message(sprintf("[brms] Full formula: %s", deparse(full_formula)))
  
  # Set model family
  family_name <- scenario$family %||% 'gaussian'
  fam <- switch(tolower(family_name),
                gaussian = gaussian(),
                student_t = brms::student(),
                {
                  message(sprintf("[brms] Unknown family '%s', defaulting to gaussian", family_name))
                  gaussian()
                })
  
  message(sprintf("[brms] Using family: %s", family_name))

  
  # === MODEL FITTING CONFIGURATION ===
  
  # Check for random effects in formula
  has_random_effects <- grepl("\\|", formula_rhs, fixed = FALSE)
  message(sprintf("[brms] Formula includes random effects: %s", has_random_effects))
  
  # Use default priors for now to ensure stability
  # prior_key <- scenario$prior_key %||% "default"
  # pri <- build_priors(prior_key, brms_config$priors, df_sub, has_re = has_random_effects)
  pri <- NULL  # Use brms default priors
  
  # Get fitting settings from configuration
  settings <- brms_config$settings
  control_list <- list(
    adapt_delta = settings$adapt_delta %||% 0.95,
    max_treedepth = settings$max_treedepth %||% 12
  )
  
  # Set backend
  backend <- settings$backend %||% 'rstan'
  if (backend == 'cmdstanr') {
    options(brms.backend = 'cmdstanr')
    message("[brms] Using cmdstanr backend")
  } else {
    message("[brms] Using rstan backend")
  }
  
  # === MODEL FITTING WITH FALLBACK STRATEGY ===
  message("[brms] Starting model fitting...")
  
  # Helper function for safe model fitting
  fit_model_safely <- function(formula_obj, prior_obj, data_obj, description) {
    message(sprintf("[brms] Attempting fit: %s", description))
    
    fit_result <- try({
      brms::brm(
        formula = formula_obj,
        data = data_obj,
        family = fam,
        prior = prior_obj,
        chains = as.integer(settings$chains %||% 2),
        iter = settings$iter %||% 1000,
        warmup = settings$warmup %||% 500,
        cores = settings$cores %||% 2,
        control = control_list,
        init = 0,
        seed = settings$seed %||% 12345,
        refresh = 100  # Show progress every 100 iterations
      )
    }, silent = FALSE)
    
    if (inherits(fit_result, "try-error")) {
      message(sprintf("[brms] Fit failed: %s", description))
      return(NULL)
    }
    
    # Check if we got valid draws
    n_draws <- tryCatch(posterior::ndraws(fit_result), error = function(e) 0)
    if (n_draws <= 0) {
      message(sprintf("[brms] Fit produced no valid draws: %s", description))
      return(NULL)
    }
    
    message(sprintf("[brms] Fit successful: %s (%d draws)", description, n_draws))
    return(fit_result)
  }
  
  # Try primary model specification
  fit <- fit_model_safely(full_formula, pri, df_sub, "Primary model")
  fallback_info <- ""
  
  # Fallback strategy if primary model fails
  if (is.null(fit) && has_random_effects) {
    message("[brms] Primary model failed, trying without random effects...")
    
    # Remove random effects from formula
    rhs_no_re <- gsub("\\+\\s*\\(1\\s*\\|[^)]+\\)", "", formula_rhs)
    rhs_no_re <- gsub("\\+\\s*\\(1\\s*\\+[^)]+\\|[^)]+\\)", "", rhs_no_re)  # Also remove random slopes
    formula_no_re <- stats::as.formula(paste(lhs, "~", rhs_no_re))
    
    # Build priors without random effects
    # pri_no_re <- build_priors(prior_key, brms_config$priors, df_sub, has_re = FALSE)
    pri_no_re <- NULL  # Use brms default priors
    
    fit <- fit_model_safely(formula_no_re, pri_no_re, df_sub, "No random effects")
    if (!is.null(fit)) {
      fallback_info <- "_noRE"
    }
  }
  
  # Final fallback: midpoint response if interval censoring fails
  if (is.null(fit) && !use_midpoint) {
    message("[brms] Trying midpoint response as final fallback...")
    
    df_sub$response_mid <- (df_sub$AAMR_lo2 + df_sub$AAMR_hi2) / 2
    lhs_mid <- 'response_mid'
    
    # Use simplified formula without random effects
    rhs_simple <- gsub("\\+\\s*\\([^)]+\\)", "", formula_rhs)
    formula_mid <- stats::as.formula(paste(lhs_mid, "~", rhs_simple))
    
    # pri_simple <- build_priors(prior_key, brms_config$priors, df_sub, has_re = FALSE)
    pri_simple <- NULL  # Use brms default priors
    
    fit <- fit_model_safely(formula_mid, pri_simple, df_sub, "Midpoint response")
    if (!is.null(fit)) {
      fallback_info <- "_midpoint"
    }
  }
  
  # Check if any model fit succeeded
  if (is.null(fit)) {
    stop("All model fitting attempts failed. Check data and model specification.")
  }
  
  message(sprintf("[brms] Model fitting completed successfully%s", fallback_info))
  
  # === MODEL SAVING AND OUTPUT ===
  message("[brms] Saving model results...")
  
  # Generate output file prefix based on scenario and timestamp
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  file_prefix <- paste0(timestamp, "_", scenario_name)
  
  # Save fitted model object
  model_dir <- results_config$model_fits
  dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)
  rds_path <- file.path(model_dir, paste0(file_prefix, ".rds"))
  saveRDS(fit, rds_path)
  message(sprintf("[brms] Saved model to: %s", rds_path))
  
  # Validate fit object
  if (!inherits(fit, 'brmsfit')) {
    stop("Invalid model fit object")
  }
  
  # Extract model summaries
  message("[brms] Extracting model summaries...")
  
  # Fixed effects
  fixed_effects <- as.data.frame(brms::fixef(fit, robust = TRUE)) %>%
    tibble::rownames_to_column("term")
  
  # Random effects (if present)
  random_effects <- NULL
  if (has_random_effects && grepl("\\|", formula_rhs, fixed = FALSE)) {
    random_effects <- tryCatch({
      as.data.frame(brms::ranef(fit))
    }, error = function(e) {
      message("[brms] Warning: Could not extract random effects")
      NULL
    })
  }
  
  # Model summary
  model_summary <- as.data.frame(summary(fit)$fixed) %>%
    tibble::rownames_to_column("term")
  
  # Save CSV outputs
  csv_dir <- results_config$csv_outputs
  dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)
  
  readr::write_csv(fixed_effects, file.path(csv_dir, paste0(file_prefix, "_fixef.csv")))
  readr::write_csv(model_summary, file.path(csv_dir, paste0(file_prefix, "_summary.csv")))
  
  if (!is.null(random_effects)) {
    suppressWarnings({
      readr::write_csv(random_effects, file.path(csv_dir, paste0(file_prefix, "_ranef.csv")))
    })
  }
  
  # Model diagnostics
  diagnostics <- tibble::tibble(
    scenario = scenario_name,
    timestamp = timestamp,
    n_obs = nrow(df_sub),
    formula_key = formula_key,
    family = family_name,
    fallback_mode = fallback_info,
    rhat_max = suppressWarnings(max(posterior::rhat(fit), na.rm = TRUE)),
    ess_bulk_min = suppressWarnings(min(posterior::ess_bulk(fit), na.rm = TRUE)),
    ess_tail_min = suppressWarnings(min(posterior::ess_tail(fit), na.rm = TRUE))
  )
  
  readr::write_csv(diagnostics, file.path(csv_dir, paste0(file_prefix, "_diagnostics.csv")))
  message(sprintf("[brms] Saved summaries to: %s", csv_dir))

  message(sprintf("[brms] Scenario '%s' completed successfully!", scenario_name))
  message("[brms] Note: Use 04_process_results.R for LMM-compatible format conversion")
}

if (identical(environment(), globalenv())) {
  main()
}
