#!/usr/bin/env Rscript
#' 
#' BRMS Model Runner - Refactored for Intelligent Self-Configuration
#' 
#' 职责反转后的设计：
#' 1. 接收单一参数：scenario_name (通过 --scenario 标志)
#' 2. 自我配置：从 config.yaml 读取场景的完整定义
#' 3. 数据准备：内置所有必要的筛选和变量选择逻辑
#' 4. 智能输出：基于场景配置生成适当的输出文件名
#' 
#' 用法: Rscript 02_run_brms_model.R --scenario LungCancer_TotalEQI_Lag5_AllRUCC
#'

suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(readr)
  library(brms)
  library(posterior)
  library(tibble)
  library(stringr)
})

# Utility operator
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Find project root by looking for config.yaml
find_project_root <- function(start = getwd()) {
  cur <- normalizePath(start)
  for (i in 1:6) {
    cand <- file.path(cur, "config.yaml")
    if (file.exists(cand)) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  stop("config.yaml not found from starting directory: ", start)
}

#' Read configuration from project root
read_config <- function() {
  # Try to find root from script location first, then from working directory
  script_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- script_args[grep("^--file=", script_args)]
  
  if (length(file_arg) > 0) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    root <- tryCatch(find_project_root(script_dir), 
                     error = function(e) find_project_root(getwd()))
  } else {
    root <- find_project_root(getwd())
  }
  
  config_path <- file.path(root, "config.yaml")
  message(sprintf("[brms] Reading config from: %s", config_path))
  yaml::read_yaml(config_path)
}

#' Parse command line arguments for scenario name
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  # Look for --scenario flag
  scenario_idx <- which(args == "--scenario")
  if (length(scenario_idx) == 0) {
    stop("Usage: Rscript 02_run_brms_model.R --scenario <scenario_name> [--config <config_file>]")
  }
  
  if (scenario_idx == length(args)) {
    stop("--scenario flag requires a scenario name")
  }
  
  scenario_name <- args[scenario_idx + 1]
  message(sprintf("[brms] Target scenario: %s", scenario_name))
  
  # Look for optional --config flag
  config_file <- NULL
  config_idx <- which(args == "--config")
  if (length(config_idx) > 0) {
    if (config_idx < length(args)) {
      config_file <- args[config_idx + 1]
      message(sprintf("[brms] Using config file: %s", config_file))
    }
  }
  
  return(list(scenario_name = scenario_name, config_file = config_file))
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

#' Resolve EQI column name based on domain and RUCC filter
resolve_eqi_column <- function(domain, rucc_filter, field_map) {
  # Domain-specific column mapping
  domain_map <- list(
    total = field_map$eqi_quintile %||% 'EQI',
    air = 'EQI_air',
    water = 'EQI_water', 
    land = 'EQI_land',
    built = 'EQI_built',
    sociodemographic = 'EQI_Sociodemographic'
  )
  
  # RUCC-stratified column mapping  
  rucc_map <- list(
    total = 'RUCC_EQI',
    air = 'RUCC_EQI_air',
    water = 'RUCC_EQI_water',
    land = 'RUCC_EQI_land', 
    built = 'RUCC_EQI_built',
    sociodemographic = 'RUCC_EQI_Sociodemographic'
  )
  
  # Normalize domain key
  domain_key <- if (is.null(domain) || domain == '' || is.na(domain)) 'total' else as.character(domain)
  if (!domain_key %in% names(domain_map)) {
    message(sprintf("[brms] Unknown domain '%s', defaulting to 'total'", domain_key))
    domain_key <- 'total'
  }
  
  # Choose column based on RUCC filter
  has_rucc_filter <- !is.null(rucc_filter) && length(rucc_filter) > 0 && !all(is.na(rucc_filter))
  
  if (has_rucc_filter) {
    column_name <- rucc_map[[domain_key]]
  } else {
    column_name <- domain_map[[domain_key]]
  }
  
  message(sprintf("[brms] Using EQI column: %s (domain=%s, rucc_filter=%s)", 
                  column_name, domain_key, ifelse(is.null(rucc_filter), "none", toString(rucc_filter))))
  
  return(column_name)
}

#' Build priors for brms model
build_priors <- function(config, data, has_re = FALSE) {
  prior_config <- config$brms_analysis$priors
  
  # Build prior list
  priors <- c(
    # Intercept prior with auto-centering
    if (prior_config$intercept$strategy == "auto_center") {
      center_val <- tryCatch({
        median(c(data$AAMR_lo2, data$AAMR_hi2), na.rm = TRUE)
      }, error = function(e) {
        prior_config$intercept$fallback_mean %||% 100.0
      })
      brms::prior(sprintf("normal(%s, %s)", 
                         center_val, 
                         prior_config$intercept$fallback_sd %||% 50.0), 
                 class = "Intercept")
    } else {
      brms::prior(sprintf("normal(%s, %s)",
                         prior_config$intercept$fallback_mean %||% 100.0,
                         prior_config$intercept$fallback_sd %||% 50.0),
                 class = "Intercept")
    },
    
    # Fixed effects priors
    brms::prior(sprintf("normal(0, %s)", prior_config$fixed_effects$class_b_sd %||% 5.0), 
               class = "b"),
    
    # Residual standard deviation
    brms::prior(sprintf("student_t(%s, 0, %s)", 
                       prior_config$sigma$df %||% 3,
                       prior_config$sigma$scale %||% 30.0), 
               class = "sigma")
  )
  
  # Add random effects priors if needed
  if (has_re) {
    re_config <- prior_config$random_effects$state_sd
    priors <- c(priors,
      brms::prior(sprintf("student_t(%s, 0, %s)",
                         re_config$df %||% 3, 
                         re_config$scale %||% 20.0),
                 class = "sd")
    )
  }
  
  return(priors)
}

#' Build model formulas with dynamic EQI column
build_formulas <- function(field_map, eqi_col_active) {
  eqi_q <- eqi_col_active %||% (field_map$eqi_quintile %||% "EQI") 
  eqi_total <- field_map$eqi_total %||% "EQI"
  state <- field_map$state %||% "State"

  # Return RHS strings to be combined into full formula later
  list(
    total_eqi_quintile = sprintf("%s + Smoking_Rate_std + (1 | %s)", eqi_q, state),
    total_eqi_quintile_nore = sprintf("%s + Smoking_Rate_std", eqi_q),
    random_slope_eqi = sprintf("scale(%s) + Smoking_Rate_std + (1 + scale(%s) | %s)", eqi_total, eqi_total, state)
  )
}

build_priors <- function(cfg, df_sub, has_re = TRUE) {
  pr <- cfg$brms_analysis$priors
  # Intercept prior center: auto_center on mid-point of observed intervals' mean
  if (!is.null(pr$intercept$strategy) && pr$intercept$strategy == "auto_center") {
    lo <- df_sub[[cfg$brms_analysis$field_map$aamr_lower %||% 'AAMR_lower']]
    hi <- df_sub[[cfg$brms_analysis$field_map$aamr_upper %||% 'AAMR_upper']]
    mid <- mean((lo + hi) / 2, na.rm = TRUE)
    icpt_mean <- if (is.finite(mid)) mid else (pr$intercept$fallback_mean %||% 100.0)
  } else {
    icpt_mean <- pr$intercept$fallback_mean %||% 100.0
  }
  icpt_sd <- pr$intercept$fallback_sd %||% 50.0

  # Fixed effects prior sd
  b_sd <- pr$fixed_effects$class_b_sd %||% 5.0

  # Sigma prior
  sig <- pr$sigma
  if (!is.null(sig) && (sig$type %||% 'student_t') == 'student_t') {
  sigma_prior <- brms::set_prior(sprintf("student_t(%s, 0, %s)", sig$df %||% 3, sig$scale %||% 30), class = "sigma")
  } else {
  sigma_prior <- brms::set_prior("normal(0, 30)", class = "sigma")
  }

  pri_list <- list(
    brms::set_prior(sprintf("normal(%s, %s)", round(icpt_mean, 2), round(icpt_sd, 2)), class = "Intercept"),
    brms::set_prior(sprintf("normal(0, %s)", b_sd), class = "b"),
    sigma_prior
  )
  if (isTRUE(has_re)) {
    re <- pr$random_effects$state_sd
    if (!is.null(re) && (re$type %||% 'student_t') == 'student_t') {
      pri_list <- c(pri_list, list(brms::set_prior(sprintf("student_t(%s, 0, %s)", re$df %||% 3, re$scale %||% 20), class = "sd")))
    } else {
      pri_list <- c(pri_list, list(brms::set_prior("normal(0, 20)", class = "sd")))
    }
  }
  do.call(c, pri_list)
}

main <- function() {
  message("[brms] Starting intelligent brms model runner (refactored)")
  
  # Parse command line arguments
  args_list <- parse_args()
  scenario_name <- args_list$scenario_name
  config_file <- args_list$config_file
  
  # Load configuration (use custom config if provided)
  if (!is.null(config_file)) {
    message(sprintf("[brms] Reading config from: %s", config_file))
    cfg <- yaml::read_yaml(config_file)
  } else {
    cfg <- read_config()
  }
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
  eqi_active_col <- resolve_eqi_column(scenario$domain, scenario$rucc_filter, field_map)

  
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
  
  # Build formulas based on scenario configuration
  formulas <- build_formulas(field_map, eqi_active_col)
  formula_key <- scenario$formula_key
  
  if (is.null(formulas[[formula_key]])) {
    available_formulas <- paste(names(formulas), collapse = ", ")
    stop(sprintf("Unknown formula_key '%s'. Available: %s", formula_key, available_formulas))
  }
  
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
  rhs <- paste("~", formulas[[formula_key]])
  full_formula <- stats::as.formula(paste(lhs, rhs))
  
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
  has_random_effects <- grepl("\\|", formulas[[formula_key]], fixed = FALSE)
  message(sprintf("[brms] Formula includes random effects: %s", has_random_effects))
  
  # Build priors
  pri <- build_priors(cfg, df_sub, has_re = has_random_effects)
  
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
    rhs_no_re <- gsub("\\+\\s*\\(1\\s*\\|[^)]+\\)", "", formulas[[formula_key]])
    rhs_no_re <- gsub("\\+\\s*\\(1\\s*\\+[^)]+\\|[^)]+\\)", "", rhs_no_re)  # Also remove random slopes
    formula_no_re <- stats::as.formula(paste(lhs, "~", rhs_no_re))
    
    # Build priors without random effects
    pri_no_re <- build_priors(cfg, df_sub, has_re = FALSE)
    
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
    rhs_simple <- gsub("\\+\\s*\\([^)]+\\)", "", formulas[[formula_key]])
    formula_mid <- stats::as.formula(paste(lhs_mid, "~", rhs_simple))
    
    pri_simple <- build_priors(cfg, df_sub, has_re = FALSE)
    
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
  if (has_random_effects && grepl("\\|", formulas[[formula_key]], fixed = FALSE)) {
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

    
  # === LMM-COMPATIBLE OUTPUT FORMAT ===
  message("[brms] Creating LMM-compatible summary output...")
  
  # Extract posterior draws and fixed effects for coefficient formatting
  draws <- posterior::as_draws_df(fit)
  fixed_effects_table <- as.data.frame(brms::fixef(fit, robust = TRUE))
  fixed_effects_table$term <- rownames(fixed_effects_table)
  
  # Identify EQI quintile coefficient names (Q2-Q5, Q1 is reference)
  find_quintile_term <- function(q) {
    possible_names <- c(
      sprintf("%sQ%s", eqi_active_col, q),
      sprintf("%s_Q%s", eqi_active_col, q),
      sprintf("EQI_quintileQ%s", q),
      sprintf("EQIQ%s", q)
    )
    
    matches <- intersect(possible_names, fixed_effects_table$term)
    if (length(matches) > 0) return(matches[1])
    return(NA_character_)
  }
  
  quintile_terms <- sapply(2:5, find_quintile_term)
  
  # Calculate significance stars based on posterior probability
  get_significance_stars <- function(term_name) {
    if (is.na(term_name)) return("")
    
    draw_column <- paste0("b_", term_name)
    if (!draw_column %in% names(draws)) return("")
    
    posterior_samples <- draws[[draw_column]]
    prob_positive <- mean(posterior_samples > 0, na.rm = TRUE)
    prob_negative <- mean(posterior_samples < 0, na.rm = TRUE)
    max_prob <- max(prob_positive, prob_negative)
    
    if (max_prob >= 0.999) return("***")
    if (max_prob >= 0.99) return("**")
    if (max_prob >= 0.95) return("*")
    return("")
  }
  
  # Format coefficient display: estimate(lower, upper)stars
  format_coefficient <- function(term_name) {
    if (is.na(term_name)) return("")
    
    term_row <- fixed_effects_table[fixed_effects_table$term == term_name, , drop = FALSE]
    if (nrow(term_row) == 0) return("")
    
    estimate <- sprintf("%.2f", term_row$Estimate)
    lower <- sprintf("%.2f", term_row$`Q2.5`)
    upper <- sprintf("%.2f", term_row$`Q97.5`)
    stars <- get_significance_stars(term_name)
    
    return(paste0(estimate, "(", lower, ", ", upper, ")", stars))
  }
  
  # Build model label using scenario configuration
  build_model_label <- function(scenario, fallback_info) {
    base_label <- switch(scenario$domain,
                        total = "EQI",
                        air = "EQI_air", 
                        water = "EQI_water",
                        land = "EQI_land",
                        built = "EQI_built",
                        sociodemographic = "EQI_Sociodemographic",
                        "EQI")  # default
    
    # Add RUCC prefix if filtering by RUCC
    if (!is.null(scenario$rucc_filter)) {
      rucc_str <- if (length(scenario$rucc_filter) == 1) {
        as.character(scenario$rucc_filter)
      } else {
        paste(range(scenario$rucc_filter), collapse = "_")
      }
      base_label <- paste0("RUCC", rucc_str, "_", base_label)
    }
    
    # Add fallback suffix
    label_suffix <- scenario$output_suffix %||% ""
    if (fallback_info != "") {
      label_suffix <- paste0(label_suffix, fallback_info)
    }
    
    return(paste0(base_label, label_suffix))
  }
  
  # Map EQI period codes for output compatibility
  format_eqi_period <- function(period_value) {
    period_str <- as.character(period_value)
    case_when(
      period_str %in% c("0005", "00_05", "2000-2005") ~ "2000_2005",
      period_str %in% c("0610", "06_10", "2006-2010") ~ "2006_2010", 
      TRUE ~ gsub("-", "_", period_str)
    )
  }
  
  # Extract key values from scenario and processed data
  cancer_code <- scenario$cancer_type
  eqi_period_formatted <- format_eqi_period(scenario$eqi_period %||% "")
  aamr_period_formatted <- gsub("-", "_", as.character(scenario$time_period %||% ""))
  lag_value <- scenario$lag_years %||% NA
  model_label <- build_model_label(scenario, fallback_info)
  
  # Create output row in LMM format
  summary_row <- tibble::tibble(
    ICD_Code = cancer_code,
    EQI_Period = eqi_period_formatted,
    AAMR_Period = aamr_period_formatted,
    Lag = lag_value,
    Model = model_label,
    Q1 = "0.00",  # Reference category
    Q2 = format_coefficient(quintile_terms[1]),
    Q3 = format_coefficient(quintile_terms[2]),
    Q4 = format_coefficient(quintile_terms[3]),
    Q5 = format_coefficient(quintile_terms[4])
  )
  
  # Write to brms-specific results file
  output_dir <- results_config$base_dir
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  output_file <- file.path(output_dir, paste0("brms_", cancer_code, ".csv"))
  
  # Append to existing file or create new one
  if (file.exists(output_file)) {
    suppressWarnings({
      readr::write_csv(summary_row, output_file, append = TRUE, col_names = FALSE)
    })
    message(sprintf("[brms] Appended results to: %s", output_file))
  } else {
    readr::write_csv(summary_row, output_file)
    message(sprintf("[brms] Created results file: %s", output_file))
  }
  
  message(sprintf("[brms] Scenario '%s' completed successfully!", scenario_name))
}

if (identical(environment(), globalenv())) {
  main()
}
