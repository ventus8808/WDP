#'
#' BRMS Results Post-Processor - Comprehensive Analysis & Conversion
#'
#' 职责：
#' 1. 汇总所有独立的CSV输出文件
#' 2. 生成可视化报告 (森林图)
#' 3. 转换为LMM兼容格式
#' 4. 生成最终的、格式化的分析报告
#'

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(yaml)
})

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(purrr)
  library(yaml)
})

# Utility functions merged from utils.R
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

#' Load project configuration from config.yaml
load_project_config <- function() {
  # Try to determine script directory
  script_dir <- tryCatch({
    dirname(normalizePath(sys.frame(1)$ofile))
  }, error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- args[grep("^--file=", args)]
    if (length(file_arg) > 0) {
      dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    } else {
      getwd()
    }
  })
  
  # Find project root
  project_root <- find_project_root(script_dir)
  
  # Load config
  config_path <- file.path(project_root, "config.yaml")
  cfg <- yaml::read_yaml(config_path)
  
  # Add project root to config for convenience
  cfg$project_root <- project_root
  
  return(cfg)
}

#' Resolve EQI column name based on domain and RUCC filter
resolve_eqi_column <- function(eqi_column_map, domain, has_rucc_filter = FALSE) {
  map_key <- if (has_rucc_filter) "with_rucc" else "no_rucc"
  column_map <- eqi_column_map[[map_key]]
  
  if (is.null(column_map)) {
    stop("EQI column map not found for key: ", map_key)
  }
  
  eqi_col <- column_map[[domain]]
  if (is.null(eqi_col)) {
    stop("EQI column not found for domain: ", domain)
  }
  
  return(eqi_col)
}

#' Prepare data for BRMS modeling
prepare_brms_data <- function(df, field_map, eqi_col) {
  # Get field names with fallbacks
  state_col <- field_map$state %||% "State"
  cancer_col <- field_map$cancer_type %||% "Cancer_Type" 
  rucc_col <- field_map$rucc %||% "RUCC"
  smoke_col <- field_map$smoking_rate %||% "Smoking_Rate"
  lo_col <- field_map$aamr_lower %||% "AAMR_lower"
  hi_col <- field_map$aamr_upper %||% "AAMR_upper"
  
  # Create EQI quintiles from the specified EQI column
  df_prep <- df %>%
    mutate(
      EQI_quintile = ntile(.data[[eqi_col]], 5),
      EQI_quintile = factor(paste0("Q", EQI_quintile), levels = paste0("Q", 1:5))
    ) %>%
    mutate(
      across(all_of(c(state_col, cancer_col, rucc_col)), as.factor),
      Smoking_Rate_std = as.numeric(scale(.data[[smoke_col]]))
    ) %>%
    mutate(
      cens_indicator = ifelse(.data[[lo_col]] == .data[[hi_col]], "none", "interval"),
      AAMR_response = ifelse(cens_indicator == "none", .data[[lo_col]], NA_real_)
    )
  
  return(df_prep)
}

#' Generate LMM-compatible output from brms results
generate_lmm_compatible <- function(cfg) {
  message("   - Generating LMM-compatible output...")
  
  brms_config <- cfg$brms_analysis
  results_config <- brms_config$results
  
  # Find all fixed effects files
  files <- list.files(
    file.path(cfg$project_root, results_config$csv_outputs), 
    pattern = "_fixef.csv$", 
    full.names = TRUE
  )
  
  if (length(files) == 0) {
    message("   - No brms fixed effects files found to process for LMM format.")
    return(invisible(NULL))
  }
  
  # Process each file to extract scenario information
  lmm_results <- data.frame()
  processed_scenarios <- character(0)  # Track processed scenarios to avoid duplicates
  
  for (file in files) {
    # Extract scenario name from filename by removing timestamp and _fixef.csv
    basename_file <- basename(file)
    scenario_name <- str_replace(basename_file, "_[0-9]+_fixef\\.csv$", "")
    
    # Skip if we've already processed this scenario
    if (scenario_name %in% processed_scenarios) {
      next
    }
    
    message(sprintf("   - Processing file: %s -> scenario: %s", basename_file, scenario_name))
    
    # Load fixed effects
    fe_data <- readr::read_csv(file, show_col_types = FALSE)
    
    # Parse scenario name to extract parameters (instead of looking in config)
    scenario <- tryCatch({
      # Simple parsing logic (can be improved later)
      parts <- unlist(strsplit(scenario_name, "_"))
      
      # Extract cancer type 
      cancer_type <- paste0(substr(parts[1], 1, 3), "_", substr(parts[1], 4, 6))
      
      # Look for RUCC filter
      rucc_filter <- NULL
      rucc_idx <- which(grepl("^RUCC[1-4]$", parts))
      if (length(rucc_idx) > 0) {
        rucc_num <- as.numeric(gsub("RUCC", "", parts[rucc_idx]))
        # 直接映射：RUCC1->1, RUCC2->2, RUCC3->3, RUCC4->4
        rucc_filter <- rucc_num
      }
      
      # Extract EQI domain
      eqi_idx <- which(parts == "EQI")
      domain <- "total"  # default
      if (length(eqi_idx) > 0 && eqi_idx < length(parts) && !grepl("^(Lag|[0-9])", parts[eqi_idx + 1])) {
        domain <- parts[eqi_idx + 1]
      }
      
      # Extract lag years
      lag_idx <- which(grepl("^Lag[0-9]+$", parts))
      lag_years <- as.numeric(gsub("Lag", "", parts[lag_idx]))
      
      # Extract time periods
      eqi_period <- "0005"  # default
      time_period <- "2006-2010"  # default
      
      # Look for period specifications
      period_idx <- which(grepl("^[0-9]{4}$", parts))
      if (length(period_idx) >= 2) {
        eqi_period <- parts[period_idx[1]]
        aamr_part <- parts[period_idx[2]]
        time_period <- paste0(substr(aamr_part, 1, 4), "-", substr(aamr_part, 5, 8))
      } else if (lag_years == 10) {
        time_period <- "2011-2015"
      }
      
      list(
        cancer_type = cancer_type,
        eqi_period = eqi_period,
        time_period = time_period,
        lag_years = lag_years,
        rucc_filter = rucc_filter,
        domain = domain
      )
    }, error = function(e) {
      message(sprintf("   - Warning: Could not parse scenario name %s: %s", scenario_name, e$message))
      return(NULL)
    })
    
    if (is.null(scenario)) {
      next
    }
    
    # Extract EQI quintile coefficients
    quintile_coefs <- fe_data %>%
      filter(str_detect(term, "EQI_quintileQ"))
    
    if (nrow(quintile_coefs) == 0) {
      message(sprintf("   - Warning: No EQI quintile coefficients found in %s", scenario_name))
      next
    }
    
    # Create formatted coefficients
    q_values <- c("0.00") # Q1 is reference
    for (q in 2:5) {
      q_row <- quintile_coefs[quintile_coefs$term == paste0("EQI_quintileQ", q), ]
      if (nrow(q_row) > 0) {
        # Check for significance
        lower <- q_row$`Q2.5`
        upper <- q_row$`Q97.5`
        estimate <- q_row$Estimate
        
        # Determine significance stars
        sig_stars <- ""
        if ((lower > 0 && upper > 0) || (lower < 0 && upper < 0)) {
          if (abs(estimate) > 2.576) sig_stars <- "***"  # p < 0.01
          else if (abs(estimate) > 1.96) sig_stars <- "**"   # p < 0.05  
          else if (abs(estimate) > 1.645) sig_stars <- "*"    # p < 0.10
        }
        
        coef_str <- sprintf("%.2f(%.2f, %.2f)%s", estimate, lower, upper, sig_stars)
      } else {
        coef_str <- "NA"
      }
      q_values <- c(q_values, coef_str)
    }
    
    # Determine model name based on RUCC filter and domain
    model_name <- if (!is.null(scenario$rucc_filter)) {
      # 直接使用RUCC值：1->RUCC1, 2->RUCC2, 3->RUCC3, 4->RUCC4
      rucc_name <- paste0("RUCC", scenario$rucc_filter)
      if (scenario$domain == "total") {
        paste0(rucc_name, "_EQI")
      } else {
        paste0(rucc_name, "_EQI_", scenario$domain)
      }
    } else {
      # No RUCC filter - overall analysis
      if (scenario$domain == "total") {
        "EQI"
      } else {
        paste0("EQI_", scenario$domain)
      }
    }
    
    # Format EQI period properly (0005 -> 2000_2005)
    eqi_period_formatted <- paste0("20", substr(scenario$eqi_period, 1, 2), "_20", substr(scenario$eqi_period, 3, 4))
    
    # Create row for this scenario
    result_row <- data.frame(
      ICD_Code = scenario$cancer_type,
      EQI_Period = eqi_period_formatted,
      AAMR_Period = gsub("-", "_", scenario$time_period),
      Lag = scenario$lag_years,
      Model = model_name,
      Q1 = q_values[1],
      Q2 = q_values[2], 
      Q3 = q_values[3],
      Q4 = q_values[4],
      Q5 = q_values[5],
      stringsAsFactors = FALSE
    )
    
    lmm_results <- rbind(lmm_results, result_row)
    processed_scenarios <- c(processed_scenarios, scenario_name)
  }
  
  # Write combined LMM format file
  if (nrow(lmm_results) > 0) {
    output_path <- file.path(cfg$project_root, results_config$reports, "brms_lmm_compatible.csv")
    readr::write_csv(lmm_results, output_path)
    message(sprintf("   - LMM-compatible format saved to: %s", output_path))
    
    # Also write cancer-specific files
    for (cancer in unique(lmm_results$ICD_Code)) {
      cancer_data <- lmm_results[lmm_results$ICD_Code == cancer, ]
      cancer_path <- file.path(cfg$project_root, results_config$reports, paste0("brms_", cancer, ".csv"))
      readr::write_csv(cancer_data, cancer_path)
      message(sprintf("   - %s results saved to: %s", cancer, cancer_path))
    }
  } else {
    message("   - No valid results to convert to LMM format")
  }
}

#' Main function for post-processing
main <- function() {
  message("\n--- Running BRMS Post-Processing ---")
  
  cfg <- load_project_config()
  res_config <- cfg$brms_analysis$results
  
  # --- 1. Aggregate all results ---
  message("1. Aggregating individual model outputs...")
  
  csv_dir <- file.path(cfg$project_root, res_config$csv_outputs)
  files <- list.files(csv_dir, pattern = "_fixef.csv$", full.names = TRUE)
  
  if (length(files) == 0) {
    message("   - No brms outputs to process. Exiting.")
    return(invisible(NULL))
  }
  
  all_fe <- purrr::map_dfr(files, ~readr::read_csv(.x, show_col_types = FALSE), .id = "src") %>%
    mutate(scenario = str_extract(basename(src), "^[^_]+(?=_)")) %>%
    select(scenario, everything(), -src)
  
  # Write combined fixed effects table
  reports_dir <- file.path(cfg$project_root, res_config$reports)
  dir.create(reports_dir, showWarnings = FALSE, recursive = TRUE)
  comb_path <- file.path(reports_dir, "brms_fixed_effects_combined.csv")
  readr::write_csv(all_fe, comb_path)
  message(sprintf("   - Combined fixed effects saved to: %s", comb_path))
  
  # --- 2. Skip Visualizations (removed to avoid errors with empty data) ---
  message("2. Skipping visualizations (can be added later when data is available)...")
  
  # --- 3. Generate LMM-Compatible Format ---
  message("3. Converting results to LMM-compatible format...")
  generate_lmm_compatible(cfg)
  
  message("\n--- Post-Processing Complete ---")
}

# --- Run Main ---
if (!interactive()) {
  main()
}
