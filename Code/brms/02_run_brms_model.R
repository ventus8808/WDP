#'
#' BRMS Model Runner - Refactored for Intelligent Self-Configuration
#'
#' 职责反转后的设计原则：
#' 1. 接收单一参数：scenario_name (通过 --scenario 标志)
#' 2. 自我配置：从 config.yaml 读取场景的完整定义
#' 3. 数据准备：内置所有必要的筛选和变量选择逻辑
#' 4. 智能输出：基于场景配置生成独立的、原子化的输出文件
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
  library(yaml)
})

# Utility functions
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

#' Parse scenario name to extract parameters
parse_scenario_name <- function(scenario_name) {
  # Example scenario names:
  # C00C97_EQI_Lag5
  # C00C97_RUCC1_EQI_air_Lag5  
  # C00C97_EQI_0610_20112015_Lag5
  
  parts <- unlist(strsplit(scenario_name, "_"))
  
  # Extract cancer type (first part)
  cancer_type <- paste0(substr(parts[1], 1, 3), "_", substr(parts[1], 4, 6))
  
  # Look for RUCC filter
  rucc_filter <- NULL
  rucc_idx <- which(grepl("^RUCC[1-4]$", parts))
  if (length(rucc_idx) > 0) {
    rucc_num <- as.numeric(gsub("RUCC", "", parts[rucc_idx]))
    # 直接使用RUCC值：RUCC1->1, RUCC2->2, RUCC3->3, RUCC4->4
    rucc_filter <- rucc_num
  }
  
  # Extract model type (EQI or EQI_domains)
  eqi_idx <- which(parts == "EQI")
  
  if (length(eqi_idx) > 0 && eqi_idx < length(parts) && parts[eqi_idx + 1] == "domains") {
    # EQI_domains model (联合模型)
    model_type <- "EQI_domains"
    formula_type <- "multi_eqi_domains"
  } else if (length(eqi_idx) > 0) {
    # EQI model (总模型)
    model_type <- "EQI"  
    formula_type <- "total_eqi"
  } else {
    stop("No EQI model type found in scenario name")
  }
  
  # Extract lag years
  lag_idx <- which(grepl("^Lag[0-9]+$", parts))
  if (length(lag_idx) == 0) stop("No Lag found in scenario name")
  lag_years <- as.numeric(gsub("Lag", "", parts[lag_idx]))
  
  # Extract time periods (if specified)
  eqi_period <- "0005"  # default
  time_period <- "2006-2010"  # default
  
  # Look for period specifications (format: 0610_20112015)
  period_idx <- which(grepl("^[0-9]{4}_[0-9]{8}$", parts))
  if (length(period_idx) > 0) {
    period_part <- parts[period_idx]
    eqi_part <- substr(period_part, 1, 4)
    aamr_part <- substr(period_part, 6, 13)
    
    eqi_period <- eqi_part
    time_period <- paste0(substr(aamr_part, 1, 4), "-", substr(aamr_part, 5, 8))
  } else if (lag_years == 10) {
    time_period <- "2011-2015"
  }
  
  return(list(
    cancer_type = cancer_type,
    eqi_period = eqi_period,
    time_period = time_period,
    lag_years = lag_years,
    rucc_filter = rucc_filter,
    model_type = model_type,
    formula_type = formula_type
  ))
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
prepare_brms_data <- function(df, field_map, scenario_params) {
  # Get field names with fallbacks
  state_col <- field_map$state %||% "State"
  cancer_col <- field_map$cancer_type %||% "Cancer_Type" 
  rucc_col <- field_map$rucc %||% "RUCC"
  smoke_col <- field_map$smoking_rate %||% "Smoking_Rate"
  lo_col <- field_map$aamr_lower %||% "AAMR_lower"
  hi_col <- field_map$aamr_upper %||% "AAMR_upper"
  
  # 筛选数据：按癌症类型、时间期、滞后年数筛选 - 与LMM对齐
  df_filtered <- df %>%
    filter(
      Cancer_Type == scenario_params$cancer_type,
      Time_Period == scenario_params$time_period,
      EQI_Period == scenario_params$eqi_period,
      Lag_Years == scenario_params$lag_years
    )
  
  # 如果有RUCC筛选，应用RUCC筛选 - 与LMM对齐
  if (!is.null(scenario_params$rucc_filter)) {
    df_filtered <- df_filtered %>%
      filter(RUCC == scenario_params$rucc_filter)
  }
  
  # 根据模型类型创建EQI变量
  if (scenario_params$model_type == "EQI") {
    # 模型1: 总EQI模型
    df_prep <- df_filtered %>%
      mutate(
        EQI_quintile = factor(paste0("Q", .data[["EQI"]]), levels = paste0("Q", 1:5))
      )
  } else if (scenario_params$model_type == "EQI_domains") {
    # 模型2: EQI细分域联合模型 
    df_prep <- df_filtered %>%
      mutate(
        EQI_air_quintile = factor(paste0("Q", .data[["EQI_air"]]), levels = paste0("Q", 1:5)),
        EQI_water_quintile = factor(paste0("Q", .data[["EQI_water"]]), levels = paste0("Q", 1:5)),
        EQI_land_quintile = factor(paste0("Q", .data[["EQI_land"]]), levels = paste0("Q", 1:5)),
        EQI_built_quintile = factor(paste0("Q", .data[["EQI_built"]]), levels = paste0("Q", 1:5)),
        EQI_sociodemographic_quintile = factor(paste0("Q", .data[["EQI_Sociodemographic"]]), levels = paste0("Q", 1:5))
      )
  }
  
  # 通用数据处理
  df_prep <- df_prep %>%
    mutate(
      across(all_of(c(state_col, cancer_col, rucc_col)), as.factor),
      Smoking_Rate_std = as.numeric(scale(.data[[smoke_col]]))
    ) %>%
    mutate(
      cens_indicator = ifelse(.data[[lo_col]] == .data[[hi_col]], "none", "interval"),
      AAMR_response = ifelse(cens_indicator == "none", .data[[lo_col]], NA_real_),
      AAMR_midpoint = (.data[[lo_col]] + .data[[hi_col]]) / 2
    )
  
  message(sprintf("   - 筛选后数据: %d rows (Cancer: %s, Period: %s, EQI: %s, Lag: %d%s)", 
                  nrow(df_prep), 
                  scenario_params$cancer_type,
                  scenario_params$time_period,
                  scenario_params$eqi_period,
                  scenario_params$lag_years,
                  ifelse(is.null(scenario_params$rucc_filter), "", paste0(", RUCC: ", scenario_params$rucc_filter))))
  
  return(df_prep)
}

# --- Core Functions ---

#' Build BRMS formula based on scenario configuration
build_formula <- function(scenario_params, field_map) {
  formula_type <- scenario_params$formula_type
  
  # Base formula structure from config
  cfg <- load_project_config()
  base_formulas <- cfg$brms_analysis$formulas
  formula_str <- base_formulas[[formula_type]] 
  
  if (is.null(formula_str)) {
    stop("Formula type not found in config: ", formula_type)
  }
  
  # Use midpoint of AAMR interval for simple regression (avoiding censoring complications)
  full_formula <- bf(
    paste0("AAMR_midpoint ~ ", formula_str),
    family = gaussian()
  )
  
  return(full_formula)
}

#' Build BRMS priors based on scenario configuration
build_priors <- function(scenario_config, prior_map) {
  prior_key <- scenario_config$prior_key %||% "default"
  
  # Simple default priors for now - can be expanded based on config structure
  priors <- c(
    prior("normal(0, 5)", class = "Intercept"),
    prior("normal(0, 2.5)", class = "b"),
    prior("student_t(3, 0, 2.5)", class = "sd")
  )
  
  return(priors)
}

#' Main execution logic
main <- function() {
  # --- 1. Setup & Configuration ---
  # Simple command line argument parsing
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 2 || args[1] != "--scenario") {
    stop("Usage: Rscript 02_run_brms_model.R --scenario SCENARIO_NAME", call. = FALSE)
  }
  
  scenario_name <- args[2]
  cfg <- load_project_config()
  brms_config <- cfg$brms_analysis
  
  # Parse scenario name to extract parameters
  scenario <- parse_scenario_name(scenario_name)
  
  message(sprintf("--- Starting BRMS analysis for scenario: %s ---", scenario_name))
  
  # --- 2. Data Loading & Preparation ---
  message("1. Loading and preparing data...")
  
  full_data <- read_csv(
    file.path(cfg$project_root, brms_config$data_file),
    show_col_types = FALSE
  )
  
  # Prepare data for brms (筛选、使用预计算五分位数、因子转换)
  df_model <- prepare_brms_data(full_data, brms_config$field_map, scenario)
  
  message(sprintf("   - Final model data size: %d rows", nrow(df_model)))
  
  # --- 3. Model Definition ---
  message("2. Building model formula and priors...")
  
  formula <- build_formula(scenario, brms_config$field_map)
  priors <- build_priors(scenario, brms_config$priors)
  
  # --- 4. Model Fitting with Auto-Parallel Configuration ---
  message("3. Fitting BRMS model with optimized parallel settings...")
  
  # Auto-detect optimal core count
  mcmc_settings <- brms_config$mcmc_settings %||% list()
  
  if (mcmc_settings$auto_cores %||% TRUE) {
    max_cores <- parallel::detectCores()
    cores_fraction <- mcmc_settings$cores_fraction %||% 0.8
    optimal_cores <- max(1, floor(max_cores * cores_fraction))
    message(sprintf("   - Auto-detected %d cores, using %d cores (%.0f%%)", 
                    max_cores, optimal_cores, cores_fraction * 100))
  } else {
    optimal_cores <- mcmc_settings$cores %||% 4
  }
  
  # MCMC settings with performance optimization
  chains <- mcmc_settings$chains %||% 4
  iter <- mcmc_settings$iter %||% 2000
  warmup <- mcmc_settings$warmup %||% 1000
  max_treedepth <- mcmc_settings$max_treedepth %||% 12
  adapt_delta <- mcmc_settings$adapt_delta %||% 0.95
  seed <- mcmc_settings$seed %||% 12345
  
  message(sprintf("   - MCMC配置: %d chains × %d iter (%d warmup) on %d cores", 
                  chains, iter, warmup, optimal_cores))
  
  fit <- NULL
  
  tryCatch({
    fit <- brm(
      formula = formula,
      data = df_model,
      prior = priors,
      chains = chains,
      iter = iter,
      warmup = warmup,
      cores = optimal_cores,
      seed = seed,
      backend = "rstan",
      # Performance optimizations
      control = list(
        max_treedepth = max_treedepth,
        adapt_delta = adapt_delta
      ),
      # No file saving for speed
      file = NULL,
      # Silence some output for cleaner progress
      silent = 0,
      refresh = 0
    )
  }, error = function(e) {
    message(paste("   - Full model failed:", e$message))
    message("   - ATTENTION: Model fitting failed. No results will be produced for this scenario.")
  })
  
  if (is.null(fit)) {
    message("--- Scenario failed: ", scenario_name, " ---")
    return()
  }
  
  # --- 5. Convert Results to LMM Format and Save ---
  message("4. Converting results to LMM format and saving...")
  
  # Extract fixed effects
  fixed_effects <- as.data.frame(fixef(fit)) %>% tibble::rownames_to_column("term")
  
  if (scenario$model_type == "EQI") {
    # 模型1: 提取总EQI系数
    quintile_coefs <- fixed_effects[str_detect(fixed_effects$term, "EQI_quintileQ"), ]
    
    if (nrow(quintile_coefs) == 0) {
      message(sprintf("   - Warning: No EQI quintile coefficients found in %s", scenario_name))
      return()
    }
    
    # Create formatted coefficients in LMM style
    q_values <- c("0.00") # Q1 is reference - 保持字符串格式
    for (q in 2:5) {
      q_row <- quintile_coefs[quintile_coefs$term == paste0("EQI_quintileQ", q), ]
      if (nrow(q_row) > 0) {
        # Check for significance using Bayesian credible intervals
        lower <- q_row$`Q2.5`
        upper <- q_row$`Q97.5`
        estimate <- q_row$Estimate
        
        # Determine significance stars based on credible interval
        sig_stars <- ""
        if ((lower > 0 && upper > 0) || (lower < 0 && upper < 0)) {
          if (abs(estimate) > 2.576) sig_stars <- "***"  # ~99% CI
          else if (abs(estimate) > 1.96) sig_stars <- "**"   # ~95% CI  
          else if (abs(estimate) > 1.645) sig_stars <- "*"    # ~90% CI
        }
        
        coef_str <- sprintf("%.2f(%.2f, %.2f)%s", estimate, lower, upper, sig_stars)
      } else {
        coef_str <- "NA"
      }
      q_values <- c(q_values, coef_str)
    }
  } else if (scenario$model_type == "EQI_domains") {
    # 模型2: EQI细分域联合模型 - 提取5个域的系数并输出5行结果
    eqi_domains <- c("air", "water", "land", "built", "sociodemographic")
    
    message("   - Extracting coefficients for 5 EQI domains from joint model...")
    
    # 为每个域单独输出一行结果
    for (domain in eqi_domains) {
      domain_coefs <- fixed_effects[str_detect(fixed_effects$term, paste0("EQI_", domain, "_quintileQ")), ]
      
      if (nrow(domain_coefs) > 0) {
        q_values_domain <- c("0.00") # Q1 is reference
        for (q in 2:5) {
          q_row <- domain_coefs[domain_coefs$term == paste0("EQI_", domain, "_quintileQ", q), ]
          if (nrow(q_row) > 0) {
            lower <- q_row$`Q2.5`
            upper <- q_row$`Q97.5`
            estimate <- q_row$Estimate
            
            sig_stars <- ""
            if ((lower > 0 && upper > 0) || (lower < 0 && upper < 0)) {
              if (abs(estimate) > 2.576) sig_stars <- "***"
              else if (abs(estimate) > 1.96) sig_stars <- "**"
              else if (abs(estimate) > 1.645) sig_stars <- "*"
            }
            
            coef_str <- sprintf("%.2f(%.2f, %.2f)%s", estimate, lower, upper, sig_stars)
          } else {
            coef_str <- "NA"
          }
          q_values_domain <- c(q_values_domain, coef_str)
        }
        
        # 为每个域创建单独的结果行，考虑RUCC分层
        if (!is.null(scenario$rucc_filter)) {
          model_name_domain <- paste0("RUCC", scenario$rucc_filter, "_EQI_", domain)
        } else {
          model_name_domain <- paste0("EQI_", domain)
        }
        
        result_row <- data.frame(
          ICD_Code = scenario$cancer_type,
          EQI_Period = paste0("20", substr(scenario$eqi_period, 1, 2), "_20", substr(scenario$eqi_period, 3, 4)),
          AAMR_Period = gsub("-", "_", scenario$time_period),
          Lag = scenario$lag_years,
          Model = model_name_domain,
          Q1 = q_values_domain[1],
          Q2 = q_values_domain[2],
          Q3 = q_values_domain[3], 
          Q4 = q_values_domain[4],
          Q5 = q_values_domain[5],
          stringsAsFactors = FALSE
        )
        
        # 保存结果
        result_file <- file.path(cfg$project_root, "Result", "brms", paste0("brms_", scenario$cancer_type, "_Results.csv"))
        
        if (file.exists(result_file)) {
          existing_data <- read_csv(result_file, show_col_types = FALSE)
          # 确保Q1列为字符类型，避免类型冲突
          existing_data$Q1 <- as.character(existing_data$Q1)
          result_row$Q1 <- as.character(result_row$Q1)
          combined_data <- bind_rows(existing_data, result_row)
          write_csv(combined_data, result_file)
          message(sprintf("   - Appended %s result to: %s", model_name_domain, result_file))
        } else {
          write_csv(result_row, result_file)
          message(sprintf("   - Created new results file with %s: %s", model_name_domain, result_file))
        }
      } else {
        message(sprintf("   - Warning: No coefficients found for domain %s", domain))
      }
    }
    
    message("--- Scenario successful: ", scenario_name, " ---")
    return()  # 提前返回，因为已经保存了所有域的结果
  }
  
  # Determine model name based on RUCC filter and model type
  model_name <- if (!is.null(scenario$rucc_filter)) {
    rucc_name <- paste0("RUCC", scenario$rucc_filter)
    paste0(rucc_name, "_", scenario$model_type)
  } else {
    scenario$model_type  # "EQI" or "EQI_domains"
  }
  
  # Format EQI period properly (0005 -> 2000_2005)
  eqi_period_formatted <- paste0("20", substr(scenario$eqi_period, 1, 2), "_20", substr(scenario$eqi_period, 3, 4))
  
  # Create result row
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
  
  # Append to cancer-specific CSV file
  results_config <- brms_config$results
  output_file <- file.path(cfg$project_root, results_config$output_dir, 
                          gsub("\\{cancer_type\\}", scenario$cancer_type, results_config$filename_template))
  
  # Create output directory if needed
  dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
  
  # Append to file (create header if file doesn't exist)
  if (!file.exists(output_file)) {
    readr::write_csv(result_row, output_file)
    message(sprintf("   - Created new results file: %s", output_file))
  } else {
    # Read existing data, bind with new row, and write back
    existing_data <- readr::read_csv(output_file, show_col_types = FALSE)
    combined_data <- rbind(existing_data, result_row)
    readr::write_csv(combined_data, output_file)
    message(sprintf("   - Appended result to: %s", output_file))
  }
  
  message(sprintf("--- Scenario successful: %s ---", scenario_name))
}

# --- Run Main ---
if (!interactive()) {
  main()
}
