#!/usr/bin/env Rscript
#'
#' BRMS Results Post-Processor - Comprehensive Analysis & Conversion
#'
#' 职责：
#' 1. 汇总所有独立的CSV输出文件
#' 2. 生成可视化报告
#' 3. 转换为LMM兼容格式（替代02_run_brms_model.R中的追加逻辑）
#' 4. 生成最终的、格式化的分析报告
#'

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(ggplot2)
  library(tibble)
  library(purrr)
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

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Generate LMM-compatible output from brms results
generate_lmm_compatible <- function(cfg) {
  res <- cfg$brms_analysis$results
  brms_config <- cfg$brms_analysis
  
  # Find all fixed effects files
  files <- list.files(res$csv_outputs, pattern = "_fixef\\.csv$", full.names = TRUE)
  if (length(files) == 0) {
    message("[brms] No fixed effects files found for LMM conversion")
    return(invisible(NULL))
  }
  
  message(sprintf("[brms] Converting %d results to LMM format", length(files)))
  
  # Process each file
  lmm_rows <- purrr::map_dfr(files, function(file_path) {
    # Extract scenario name from filename - handle both old and new naming schemes
    filename <- basename(file_path)
    
    # New scheme: timestamp_ScenarioName_fixef.csv
    new_match <- stringr::str_match(filename, "^\\d{8}_\\d{6}_(.*)_fixef\\.csv$")
    if (!is.na(new_match[,2])) {
      scenario_name <- new_match[,2]
    } else {
      # Old scheme: timestamp_details_fixef.csv (fallback)
      old_match <- stringr::str_match(filename, "_(.*)_fixef\\.csv$")
      scenario_name <- if (!is.na(old_match[,2])) old_match[,2] else "Unknown"
    }
    
    # Find corresponding scenario config
    scenario <- NULL
    for (sc in brms_config$scenarios) {
      if (sc$name == scenario_name) {
        scenario <- sc
        break
      }
    }
    
    if (is.null(scenario)) {
      message(sprintf("[brms] Warning: No scenario config found for %s", scenario_name))
      return(NULL)
    }
    
    # Read fixed effects
    fixef_data <- readr::read_csv(file_path, show_col_types = FALSE)
    
    # Extract EQI quintile coefficients (Q2-Q5, Q1 is reference)
    eqi_col_map <- brms_config$eqi_column_map
    has_rucc <- !is.null(scenario$rucc_filter)
    map_key <- if (has_rucc) "with_rucc" else "no_rucc"
    eqi_col <- eqi_col_map[[map_key]][[scenario$domain %||% "total"]]
    
    find_quintile_coef <- function(q) {
      possible_patterns <- c(
        paste0(eqi_col, "Q", q),
        paste0("Q", q)
      )
      
      for (pattern in possible_patterns) {
        matches <- which(grepl(pattern, fixef_data$term, fixed = TRUE))
        if (length(matches) > 0) return(fixef_data[matches[1], ])
      }
      return(NULL)
    }
    
    # Format coefficient: estimate(lower, upper)
    format_coef <- function(coef_row) {
      if (is.null(coef_row)) return("")
      
      est <- sprintf("%.2f", coef_row$Estimate)
      lo <- sprintf("%.2f", coef_row$`Q2.5`)
      hi <- sprintf("%.2f", coef_row$`Q97.5`)
      
      # Add significance stars based on credible interval
      stars <- ""
      if (sign(coef_row$`Q2.5`) == sign(coef_row$`Q97.5`)) {
        if (abs(coef_row$`Q2.5`) > 0.5 && abs(coef_row$`Q97.5`) > 0.5) stars <- "***"
        else if (abs(coef_row$`Q2.5`) > 0.2 && abs(coef_row$`Q97.5`) > 0.2) stars <- "**"
        else stars <- "*"
      }
      
      return(paste0(est, "(", lo, ", ", hi, ")", stars))
    }
    
    # Build model label
    base_label <- switch(scenario$domain %||% "total",
                        total = "EQI",
                        air = "EQI_air", 
                        water = "EQI_water",
                        land = "EQI_land",
                        built = "EQI_built",
                        sociodemographic = "EQI_Sociodemographic",
                        "EQI")
    
    if (has_rucc) {
      rucc_str <- if (length(scenario$rucc_filter) == 1) {
        as.character(scenario$rucc_filter)
      } else {
        paste(range(scenario$rucc_filter), collapse = "_")
      }
      base_label <- paste0("RUCC", rucc_str, "_", base_label)
    }
    
    model_label <- paste0(base_label, scenario$output_suffix %||% "")
    
    # Format periods
    format_period <- function(period) {
      if (is.null(period)) return("")
      gsub("-", "_", as.character(period))
    }
    
    # Create LMM-compatible row
    tibble::tibble(
      ICD_Code = scenario$cancer_type,
      EQI_Period = format_period(scenario$eqi_period),
      AAMR_Period = format_period(scenario$time_period),
      Lag = scenario$lag_years %||% NA,
      Model = model_label,
      Q1 = "0.00",  # Reference
      Q2 = format_coef(find_quintile_coef(2)),
      Q3 = format_coef(find_quintile_coef(3)),
      Q4 = format_coef(find_quintile_coef(4)),
      Q5 = format_coef(find_quintile_coef(5))
    )
  })
  
  if (is.null(lmm_rows) || nrow(lmm_rows) == 0) {
    message("[brms] No LMM-compatible rows generated")
    return(invisible(NULL))
  }
  
  # Write consolidated LMM-compatible file
  output_dir <- res$base_dir
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  output_file <- file.path(output_dir, "brms_lmm_compatible.csv")
  
  readr::write_csv(lmm_rows, output_file)
  message(sprintf("[brms] Wrote LMM-compatible results to: %s", output_file))
  
  return(lmm_rows)
}

main <- function() {
  message("[brms] Starting results post-processing...")
  
  cfg <- load_project_config()
  res <- cfg$brms_analysis$results
  out_dir <- res$csv_outputs

  files <- list.files(out_dir, pattern = "_fixef\\.csv$", full.names = TRUE)
  if (length(files) == 0) {
    message("[brms] No brms outputs to process.")
    return(invisible(NULL))
  }

  # Combine all fixed effects results
  all_fe <- purrr::map_dfr(files, readr::read_csv, show_col_types = FALSE, .id = "src")
  all_fe <- all_fe %>%
    mutate(
      scenario = {
        filename <- basename(src)
        # New scheme: timestamp_ScenarioName_fixef.csv
        new_match <- stringr::str_match(filename, "^\\d{8}_\\d{6}_(.*)_fixef\\.csv$")
        ifelse(!is.na(new_match[,2]), new_match[,2], 
               # Old scheme fallback
               stringr::str_match(filename, "_(.*)_fixef\\.csv$")[,2])
      }
    ) %>%
    select(scenario, dplyr::everything(), -src)

  # Write combined fixed effects table
  comb_path <- file.path(res$reports, "brms_fixed_effects_combined.csv")
  dir.create(res$reports, showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(all_fe, comb_path)

  # Generate forest plot for main effects (excluding Intercept)
  p <- all_fe %>%
    filter(.data$term != "Intercept") %>%
    ggplot(aes(x = .data$Estimate, y = .data$term, xmin = .data$`Q2.5`, xmax = .data$`Q97.5`)) +
    geom_point() +
    geom_errorbarh(height = 0.2) +
    facet_wrap(~ .data$scenario, scales = "free_y") +
    theme_minimal() +
    labs(title = "brms fixed effects (95% CrI)", x = "Effect", y = "Term")

  fig_path <- file.path(cfg$brms_analysis$results$figures, "brms_fixed_effects_forest.png")
  dir.create(dirname(fig_path), showWarnings = FALSE, recursive = TRUE)
  ggsave(fig_path, p, width = 12, height = 8, dpi = 150)
  
  # Generate LMM-compatible output
  lmm_results <- generate_lmm_compatible(cfg)
  
  message(sprintf("[brms] Post-processing complete!"))
  message(sprintf("[brms] - Combined table: %s", comb_path))
  message(sprintf("[brms] - Forest plot: %s", fig_path))
  if (!is.null(lmm_results)) {
    message(sprintf("[brms] - LMM-compatible: %s", file.path(res$base_dir, "brms_lmm_compatible.csv")))
  }
}

if (identical(environment(), globalenv())) {
  main()
}
