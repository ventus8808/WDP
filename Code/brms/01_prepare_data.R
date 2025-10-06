#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
  library(readr)
  library(dplyr)
})

script_dir <- function() {
  # Try to infer script directory from --file
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- args[grep("^--file=", args)]
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  # Fallback to current working directory
  getwd()
}

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

read_config <- function() {
  root <- tryCatch(find_project_root(script_dir()), error = function(e) find_project_root(getwd()))
  yaml::read_yaml(file.path(root, "config.yaml"))
}

# safely get field name from map
get_field <- function(map, key, fallback) {
  val <- map[[key]]
  if (is.null(val)) fallback else val
}

main <- function() {
  cfg <- read_config()
  br <- cfg$brms_analysis
  stopifnot(!is.null(br))

  data_file <- br$data_file
  df <- readr::read_csv(data_file, show_col_types = FALSE)

  fm <- br$field_map
  state_col <- get_field(fm, 'state', 'State')
  cancer_col <- get_field(fm, 'cancer_type', 'Cancer_Type')
  rucc_col <- get_field(fm, 'rucc', 'RUCC')
  eqi_q_col <- get_field(fm, 'eqi_quintile', 'EQI_quintile')
  smoke_col <- get_field(fm, 'smoking_rate', 'Smoking_Rate')
  lo_col <- get_field(fm, 'aamr_lower', 'AAMR_lower')
  hi_col <- get_field(fm, 'aamr_upper', 'AAMR_upper')

  # 数据已包含预计算的EQI五分位数（1-5整数）- 与LMM对齐
  # 检验EQI五分位数列是否存在并且格式正确
  eqi_columns_to_check <- c('EQI', 'EQI_air', 'EQI_water', 'EQI_land', 'EQI_built', 'EQI_Sociodemographic',
                           'RUCC_EQI', 'RUCC_EQI_air', 'RUCC_EQI_water', 'RUCC_EQI_land', 'RUCC_EQI_built', 'RUCC_EQI_Sociodemographic')
  
  for (eqi_col in eqi_columns_to_check) {
    if (eqi_col %in% names(df)) {
      # 检查非缺失值是否都在1-5范围内
      non_na_values <- df[[eqi_col]][!is.na(df[[eqi_col]])]
      if (length(non_na_values) > 0 && all(non_na_values %in% 1:5)) {
        na_count <- sum(is.na(df[[eqi_col]]))
        cat(sprintf("✓ %s: 有效五分位数 (1-5), %d个缺失值\n", eqi_col, na_count))
      } else {
        warning(sprintf("⚠ %s column contains invalid values (not in 1-5 range)", eqi_col))
      }
    }
  }

  # Factor only columns that exist
  factor_cols <- intersect(c(state_col, cancer_col, rucc_col, eqi_q_col), names(df))

  df2 <- df %>%
    mutate(
      across(all_of(factor_cols), as.factor),
      Smoking_Rate_std = as.numeric(scale(.data[[smoke_col]]))
    ) %>%
    mutate(cens_indicator = ifelse(.data[[lo_col]] == .data[[hi_col]], "none", "interval")) %>%
    mutate(AAMR_response = ifelse(.data[["cens_indicator"]] == "none", .data[[lo_col]], NA_real_))

  # Quick validity checks
  n_bad <- sum(df2[[lo_col]] > df2[[hi_col]], na.rm = TRUE)
  if (n_bad > 0) {
    warning(sprintf("Found %d rows with AAMR_lower > AAMR_upper", n_bad))
  }

  # Optional: write a snapshot for audit under Data/Processed
  proc_dir <- cfg$data_directories$processed %||% 'Data/Processed'
  out_dir <- file.path(proc_dir, 'EQI')
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_file <- file.path(out_dir, 'EQI_AAMR_Interval_AnalysisReady.csv')
  readr::write_csv(df2, out_file)
  cat(sprintf("Analysis-ready snapshot written to %s (rows=%d)\n", out_file, nrow(df2)))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

if (identical(environment(), globalenv())) {
  main()
}
