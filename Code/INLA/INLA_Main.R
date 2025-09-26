#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# WDP BYM INLA Production Analysis Script
# Production-level pesticide exposure and health outcome analysis
# Generates 24 results per category: 3 estimates × 2 lag periods × 4 models
# Author: WDP Analysis Team
# Date: 2024

# Suppress warnings and load essential libraries
options(warn = -1)
suppressMessages({
  # Load here package for robust path management
  if (!require(here, quietly = TRUE)) {
    install.packages("here")
    library(here, quietly = TRUE)
  }
  
  library(INLA, quietly = TRUE)
  library(dplyr, quietly = TRUE)
  library(readr, quietly = TRUE)
  library(yaml, quietly = TRUE)
  library(argparse, quietly = TRUE)
  library(progress, quietly = TRUE)
})

# Print project root for verification
cat(sprintf("📁 Project root identified at: %s\n", here()))

# Configure INLA environment first
Sys.setenv(INLA_DEBUG = "0")  # Reduce INLA debug output
Sys.setenv(INLA_HOME = system.file(package = "INLA"))

# Prefer node-local temporary directory (SLURM_TMPDIR or system TMPDIR), fallback to project temp
tmp_base <- if (nzchar(Sys.getenv("SLURM_TMPDIR"))) {
  Sys.getenv("SLURM_TMPDIR")
} else if (nzchar(Sys.getenv("TMPDIR"))) {
  Sys.getenv("TMPDIR")
} else {
  here("temp")
}
project_temp <- file.path(tmp_base, "inla_temp")
if (!dir.exists(project_temp)) {
  dir.create(project_temp, recursive = TRUE, showWarnings = FALSE)
  cat(sprintf("📁 Created INLA temp directory: %s\n", project_temp))
}

# Set TMP env vars to the selected directory
Sys.setenv(TMPDIR = project_temp)
Sys.setenv(TMP = project_temp)
Sys.setenv(TEMP = project_temp)
cat(sprintf("🧊 Using temporary directory: %s\n", project_temp))

# Check INLA installation
if (!require("INLA", character.only = TRUE, quietly = TRUE)) {
  stop("❌ INLA package not found. Please install INLA first.")
}

# Create temporary directories
temp_dir <- here("temp")
graphs_dir <- file.path(temp_dir, "graphs")
models_dir <- file.path(temp_dir, "models")

for (dir in c(temp_dir, graphs_dir, models_dir)) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    cat(sprintf("📁 Created directory: %s\n", dir))
  }
}

# Configure INLA settings step by step
tryCatch({
  inla.setOption(verbose = FALSE)
  # Respect SLURM CPU allocation when available
  cpus <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")))
  threads <- if (!is.na(cpus) && cpus > 0) paste0(cpus, ":1") else "4:1"
  inla.setOption(num.threads = threads)
  inla.setOption(inla.mode = "experimental")  # Use experimental mode for better stability
  # Force INLA to keep and use a local working directory on node-local storage
  local_work_dir <- file.path(project_temp, "inla_work")
  if (!dir.exists(local_work_dir)) dir.create(local_work_dir, recursive = TRUE, showWarnings = FALSE)
  inla.setOption(keep = TRUE)
  inla.setOption(working.directory = local_work_dir)
  inla.setOption(safe = TRUE)
  # Diagnostics: show INLA options and binary path
  iwdir <- inla.getOption("working.directory")
  ickeep <- inla.getOption("keep")
  ithreads <- inla.getOption("num.threads")
  icall <- inla.getOption("inla.call")
  cat(sprintf("✅ INLA basic configuration complete (threads=%s, workdir=%s, keep=%s)\n", ithreads, iwdir, as.character(ickeep)))
  cat(sprintf("🔎 INLA binary: %s (exists=%s)\n", icall, file.exists(icall)))
}, error = function(e) {
  cat(sprintf("⚠️ INLA configuration warning: %s\n", e$message))
})

cat("🔧 INLA configured for HPC environment\n")

cat("WDP BYM INLA Production Analysis System\n")
cat("======================================\n")

source("INLA_Utils/INLA_Utils_Data.R")
source("INLA_Utils/INLA_Utils_Model.R")
source("INLA_Utils/INLA_Utils_Results.R")
source("INLA_Utils/INLA_Utils_Dashboard.R")
source("INLA_Utils/INLA_Utils_Validation.R")
source("INLA_Utils/INLA_Utils_Logger.R")

# 全局debug开关
debug <- TRUE
if (exists("config")) config$logging$level <- "DEBUG"

#' Parse command line arguments
#' @return List of parsed arguments
parse_command_args <- function() {
  parser <- ArgumentParser(description = 'WDP BYM INLA Production Analysis')

  # Input parameters
  parser$add_argument('--config', type = 'character',
                      default = 'INLA_Config/analysis_config.yaml',
                      help = 'Path to configuration YAML file')
  parser$add_argument('--disease-code', type = 'character', default = 'C81-C96',
                      help = 'Disease code to analyze')
  parser$add_argument('--measure-type', type = 'character', default = 'Weight',
                      help = 'Comma-separated measure types: Weight,Density')
  parser$add_argument('--pesticide-category', type = 'character', default = 'TEST',
                      help = 'Pesticide category: specific name, ALL, TOP5, or TEST')
  parser$add_argument('--estimate-types', type = 'character', default = 'avg',
                      help = 'Comma-separated estimate types: min,avg,max or subset')
  parser$add_argument('--lag-years', type = 'character', default = '5',
                      help = 'Comma-separated lag years: 5,10 or subset')
  parser$add_argument('--model-types', type = 'character', default = 'M0,M1,M2,M3',
                      help = 'Comma-separated model types: M0,M1,M2,M3 or subset')
  parser$add_argument('--output-file', type = 'character', default = '',
                      help = 'Custom output filename (optional)')
  parser$add_argument('--verbose', action = 'store_true', default = FALSE,
                      help = 'Enable verbose logging')
  parser$add_argument('--dry-run', action = 'store_true', default = FALSE,
                      help = 'Validate inputs without running analysis')

  if (interactive()) {
    # Default values for interactive testing
    return(list(
      config = 'INLA_Config/analysis_config.yaml',
      disease_code = 'C81-C96',
      measure_type = 'Weight',
      pesticide_category = 'TEST',
      estimate_types = 'avg',
      lag_years = '5',
      model_types = 'M0,M1,M2,M3',
      output_file = '',
      verbose = TRUE,
      dry_run = FALSE
    ))
  } else {
    args <- parser$parse_args()
    # Convert hyphenated arguments to underscores
    names(args) <- gsub("-", "_", names(args))
    return(args)
  }
}

#' Create necessary directory structure
#' @param config Configuration list
create_directory_structure <- function(config) {
  cat("📁 Creating necessary directory structure...\n")
  
  # Create output directory
  output_dir <- here(config$output$base_dir)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    cat(sprintf("  ✓ Created output directory: %s\n", output_dir))
  }
  
  # Create temp directories using configuration
  temp_dirs <- c(
    here(config$model_fitting$spatial$graph_dir),
    here("Code/INLA/INLA_Temp/models"),  # 保留models和logs目录
    here("Code/INLA/INLA_Temp/logs")
  )
  
  for (temp_dir in temp_dirs) {
    if (!dir.exists(temp_dir)) {
      dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
      cat(sprintf("  ✓ Created temp directory: %s\n", basename(temp_dir)))
    }
  }
  
  cat("✓ Directory structure ready\n")
}

#' Load and validate configuration
#' @param config_path Path to YAML configuration file
#' @return Configuration list
load_config <- function(config_path) {
  if (!file.exists(config_path)) {
    stop(sprintf("Configuration file not found: %s", config_path))
  }

  config <- read_yaml(config_path)
  cat(sprintf("✓ Configuration loaded from: %s\n", config_path))

  # Validate essential configuration sections
  required_sections <- c("data_paths", "analysis", "output", "model_fitting")
  missing_sections <- setdiff(required_sections, names(config))

  if (length(missing_sections) > 0) {
    stop(sprintf("Missing configuration sections: %s",
                 paste(missing_sections, collapse = ", ")))
  }

  return(config)
}

#' Parse comma-separated argument strings
#' @param arg_string Comma-separated string
#' @param valid_options Vector of valid options (optional)
#' @return Vector of parsed values
parse_list_argument <- function(arg_string, valid_options = NULL) {
  parsed <- trimws(strsplit(arg_string, ",")[[1]])

  if (!is.null(valid_options)) {
    invalid <- setdiff(parsed, valid_options)
    if (length(invalid) > 0) {
      stop(sprintf("Invalid options: %s. Valid options: %s",
                   paste(invalid, collapse = ", "),
                   paste(valid_options, collapse = ", ")))
    }
  }

  return(parsed)
}

#' Generate output filename
#' @param config Configuration list
#' @param disease_code Disease code
#' @param exposure_name Name of the compound or category (optional)
#' @param custom_filename Custom filename from arguments (optional)
#' @return Full output file path
generate_output_filename <- function(config, disease_code, exposure_name = NULL, custom_filename = "") {
  # If a specific output file is requested, use it directly.
  if (!is.null(custom_filename) && custom_filename != "") {
    if (!grepl("\\.csv$", custom_filename)) {
      custom_filename <- paste0(custom_filename, ".csv")
    }
    return(file.path(config$output$base_dir, custom_filename))
  }

  # Otherwise, generate a filename from a template.
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

  # Sanitize exposure name for filename
  safe_exposure_name <- if (!is.null(exposure_name)) {
    gsub("[^a-zA-Z0-9_.-]", "_", exposure_name)
  } else {
    "Analysis"
  }

  # Replace placeholders in the template
  filename <- config$output$filename_template
  filename <- gsub("\\{disease\\}", disease_code, filename)
  filename <- gsub("\\{exposure\\}", safe_exposure_name, filename) # New placeholder
  filename <- gsub("\\{timestamp\\}", timestamp, filename)

  return(file.path(config$output$base_dir, filename))
}

#' Process a single combination of parameters
#' @param data_list Loaded data
#' @param analysis_info Analysis parameters
#' @param config Configuration
#' @return List with success status and result row
process_single_combination <- function(data_list, analysis_info, config) {
  start_time <- Sys.time()

  tryCatch({
    # Build pesticide column name
    pesticide_col <- sprintf("cat%d_%s",
                            analysis_info$category_id,
                            analysis_info$estimate_type)

    # Check if pesticide column exists
    if (!check_pesticide_column(data_list$pesticide, pesticide_col)) {
      return(list(
        success = FALSE,
        result = create_failed_result_row(analysis_info, "Column not found", config),
        error = sprintf("Column %s not found", pesticide_col)
      ))
    }

    # Prepare model data
    model_data <- prepare_model_data(
      data_list, pesticide_col, analysis_info$lag_years, config)
    # 保存模型数据以便debug（使用正确的路径）
    debug_path <- file.path(project_temp, "debug_model_data.rds")
    saveRDS(model_data, file = debug_path)

    # Check minimum sample size
    if (nrow(model_data) < config$data_processing$min_thresholds$records_per_analysis) {
      return(list(
        success = FALSE,
        result = create_failed_result_row(analysis_info, "Insufficient data", config),
        error = sprintf("Insufficient data: %d records", nrow(model_data))
      ))
    }

    # Create spatial structure
    cat("  📊 About to create spatial structure...\n")
    counties_in_data <- unique(model_data$COUNTY_FIPS)
    
    # <<<<<<<  新增诊断代码开始  >>>>>>>
    cat("\n--- [DEBUG] Checking counties_in_data ---\n")
    cat(sprintf("Number of unique counties found: %d\n", length(counties_in_data)))
    cat("First 10 county FIPS:\n")
    print(head(counties_in_data, 10))
    cat("Structure of the vector:\n")
    print(str(counties_in_data))
    cat("--- [DEBUG] End of check ---\n\n")
    # <<<<<<<  新增诊断代码结束  >>>>>>>

    # 添加错误处理
    if(length(counties_in_data) == 0) {
      cat("  ❌ Error: counties_in_data is empty!\n")
      return(list(
        success = FALSE,
        result = create_failed_result_row(analysis_info, "No counties in data", config),
        error = "No counties in data"
      ))
    }
    
    spatial_structure <- tryCatch({
      create_spatial_structure(
        data_list$adjacency, counties_in_data, analysis_info$category_id, config)
    }, error = function(e) {
      cat(sprintf("  ❌ Error in create_spatial_structure: %s\n", e$message))
      return(NULL)
    })
    
    if(is.null(spatial_structure)) {
      cat("  ❌ Error: spatial_structure is NULL\n")
      return(list(
        success = FALSE,
        result = create_failed_result_row(analysis_info, "Spatial structure creation failed", config),
        error = "Spatial structure creation failed"
      ))
    }
      
    cat("  📊 Spatial structure created successfully\n")

    # Build model formula
    cat("  📊 About to build model formula...\n")
    formula <- tryCatch({
      build_model_formula(
        analysis_info$model_type, model_data, spatial_structure$filepath, config)
    }, error = function(e) {
      cat(sprintf("  ❌ Error in build_model_formula: %s\n", e$message))
      return(NULL)
    })
    
    if(is.null(formula)) {
      cat("  ❌ Error: formula is NULL\n")
      return(list(
        success = FALSE,
        result = create_failed_result_row(analysis_info, "Model formula creation failed", config),
        error = "Model formula creation failed"
      ))
    }
    cat("  📊 Model formula built successfully\n")

    # Fit INLA model
    cat("  📊 About to fit INLA model...\n")
    model <- tryCatch({
      fit_inla_model(formula, model_data, config)
    }, error = function(e) {
      cat(sprintf("  ❌ Error in fit_inla_model: %s\n", e$message))
      return(NULL)
    })
    
    if(is.null(model)) {
      cat("  ❌ Error: model is NULL\n")
      return(list(
        success = FALSE,
        result = create_failed_result_row(analysis_info, "Model fitting failed", config),
        error = "Model fitting failed"
      ))
    }
    cat("  📊 INLA model fitting completed\n")

    # Create result row
    result_row <- create_result_row(analysis_info, model, model_data, config)

    # Clean up temporary files
    cleanup_spatial_files(spatial_structure$filepath)

    # Calculate processing time
    end_time <- Sys.time()
    processing_time <- as.numeric(difftime(end_time, start_time, units = "mins"))

    return(list(
      success = TRUE,
      result = result_row,
      processing_time = processing_time
    ))

  }, error = function(e) {
    end_time <- Sys.time()
    processing_time <- as.numeric(difftime(end_time, start_time, units = "mins"))

    return(list(
      success = FALSE,
      result = create_failed_result_row(analysis_info, e$message, config),
      error = e$message,
      processing_time = processing_time
    ))
  })
}

#' Analyze a single compound
#' @param compound_info Compound information (compound_id, compound_name, etc.)
#' @param args Command line arguments
#' @param config Configuration
#' @param output_path Output file path
#' @param pb Progress bar object
#' @return List of failed combinations
analyze_single_compound <- function(compound_info, args, config, output_path, pb, preloaded_data) {

  # Parse analysis parameters from arguments
  estimate_types <- parse_list_argument(args$estimate_types, config$analysis$exposure_estimates)
  lag_years <- as.numeric(parse_list_argument(args$lag_years))
  model_types <- parse_list_argument(args$model_types, names(config$analysis$models))
  measure_types <- parse_list_argument(args$measure_type, c('Weight', 'Density'))

  # Use the category ID for data column selection, but report the compound name
  category_id <- compound_info$cat_id
  failed_combinations <- list()

  # Main analysis loop: iterates through each combination of parameters
  for (measure in measure_types) {
    # Use pre-loaded data for the current measure type
    data_list <- preloaded_data[[measure]]

    if (is.null(data_list)) {
      # If data loading fails for a measure, skip all its combinations
      num_skipped <- length(estimate_types) * length(lag_years) * length(model_types)
      for (i in 1:num_skipped) {
        pb$tick(tokens = list(what = sprintf("Skipping %s: Data loading failed", measure)))
      }
      failed_combinations <- c(failed_combinations,
                               list(sprintf("All combinations for measure '%s'", measure)))
      next # Move to the next measure type
    }

    for (estimate in estimate_types) {
      for (lag in lag_years) {
        for (model in model_types) {
          # Advance the progress bar and update its label
          pb$tick(tokens = list(what = sprintf("%s/%s/%d-yr/%s", measure, estimate, lag, model)))

          # Assemble all parameters for this specific analysis run
          analysis_info <- list(
            disease_code = args$disease_code,
            measure_type = measure,
            exposure_name = compound_info$compound_name,
            category_name = compound_info$category,
            category_id = category_id,
            estimate_type = estimate,
            lag_years = lag,
            model_type = model
          )

          # Execute the model fitting and result extraction
          result <- process_single_combination(data_list, analysis_info, config)

          # Print the result of this combination to the console dashboard
          print_results_table_row(result$result)

          # Immediately write the result row to the CSV file
          write_result_row(result$result, output_path, config)

          if (!result$success) {
            # Collect information about failed runs for the final summary
            failed_combinations <- c(failed_combinations, list(result$result))
          }
        }
      }
    }
  }

  return(failed_combinations)
}

#' Analyze a single category
#' @param category_info Category information
#' @param args Command line arguments
#' @param config Configuration
#' @param output_path Output file path
#' @param preloaded_data Pre-loaded data for all measure types
#' @return List with summary statistics
analyze_single_category <- function(category_info, args, config, output_path, preloaded_data) {

  # Generate a unique output filename for this category
  output_path <- generate_output_filename(config, args$disease_code, category_info$category)

  # Initialize the output file for this specific category
  initialize_output_file(output_path, config)

  cat(sprintf("\n🔄 Processing: %s (ID: %d)\n",
              category_info$category, category_info$cat_id))
  cat(sprintf("  📄 Output will be saved to: %s\n", output_path))

  # Parse analysis parameters
  estimate_types <- parse_list_argument(args$estimate_types,
                                        config$analysis$exposure_estimates)
  lag_years <- as.numeric(parse_list_argument(args$lag_years))
  model_types <- parse_list_argument(args$model_types,
                                     names(config$analysis$models))
  # New: Get measure types from argument
  measure_types <- parse_list_argument(args$measure_type, c('Weight', 'Density'))

  # Calculate total combinations
  total_combinations <- length(measure_types) * length(estimate_types) * length(lag_years) * length(model_types)
  cat(sprintf("  📊 Running %d combinations (%d measures × %d estimates × %d lags × %d models)\n",
              total_combinations, length(measure_types), length(estimate_types), length(lag_years), length(model_types)))

  # Initialize counters
  successful <- 0
  failed <- 0
  combination_count <- 0

  # Quadruple nested loop: measure × estimates × lag years × models
  for (measure in measure_types) {
    # Use pre-loaded data for the specific measure type
    data_list <- preloaded_data[[measure]]
    if (is.null(data_list)) {
      cat(sprintf("  ❌ Skipping measure '%s' due to data loading failure.\n", measure))
      failed <- failed + (total_combinations / length(measure_types)) # Increment failure count
      next
    }

    for (estimate in estimate_types) {
      for (lag in lag_years) {
        for (model in model_types) {
          combination_count <- combination_count + 1

          # Create analysis info
          analysis_info <- list(
            disease_code = args$disease_code,
            measure_type = measure,
            exposure_name = category_info$category,
            category_name = category_info$category,
            category_id = category_info$cat_id,
            estimate_type = estimate,
            lag_years = lag,
            model_type = model
          )

          if (args$verbose) {
            cat(sprintf("    [%d/%d] %s, %s, %dy lag, %s model... ",
                        combination_count, total_combinations, measure, estimate, lag, model))
          }

          # Process single combination
          result <- process_single_combination(data_list, analysis_info, config)

          # Write result immediately
          write_result_row(result$result, output_path, config)

          if (result$success) {
            successful <- successful + 1
            if (args$verbose) {
              cat("✅\n")
            }
          } else {
            failed <- failed + 1
            if (args$verbose) {
              cat(sprintf("❌ (%s)\n", result$error))
            }
          }

          # Print detailed summary for first combination or failures
          if (combination_count == 1 || !result$success) {
            print_result_summary(result$result)
          }
        }
      }
    }
  }

  # Category summary
  success_rate <- successful / total_combinations * 100
  cat(sprintf("  📈 Category complete: %d/%d successful (%.1f%%)\n",
              successful, total_combinations, success_rate))

  return(list(
    total = total_combinations,
    successful = successful,
    failed = failed,
    success_rate = success_rate
  ))
}

# 注释：删除了过时的validate_result_row函数，
# 现在使用INLA_Utils_Results.R中的正确版本

#' Main execution function
#' @return Exit status
main <- function() {
  analysis_start_time <- Sys.time()

  cat("Starting WDP BYM INLA Production Analysis\n")
  cat("========================================\n")

  # Parse arguments and load configuration
  args <- parse_command_args()
  config <- load_config(args$config)
  
  # Override graph_dir to node-local tmp to avoid NFS issues
  local_graph_dir <- file.path(project_temp, "graphs")
  if (!dir.exists(local_graph_dir)) dir.create(local_graph_dir, recursive = TRUE, showWarnings = FALSE)
  if (!is.null(config$model_fitting) && !is.null(config$model_fitting$spatial)) {
    config$model_fitting$spatial$graph_dir <- local_graph_dir
  }
  cat(sprintf("🗂️  Spatial graphs directory set to: %s\n", local_graph_dir))
  
  # Create necessary directory structure
  create_directory_structure(config)

  # Check if we're analyzing compounds or categories
  analyze_compounds <- grepl("^compound:", args$pesticide_category)

  if (args$verbose) {
    cat("Analysis Parameters:\n")
    cat(sprintf("  Disease Code: %s\n", args$disease_code))
    cat(sprintf("  Measure Type: %s\n", args$measure_type))
    cat(sprintf("  Pesticide Category/Compound: %s\n", args$pesticide_category))
    cat(sprintf("  Estimate Types: %s\n", args$estimate_types))
    cat(sprintf("  Lag Years: %s\n", args$lag_years))
    cat(sprintf("  Model Types: %s\n", args$model_types))
    cat(sprintf("  Verbose: %s\n", args$verbose))
    cat(sprintf("  Dry Run: %s\n", args$dry_run))
    cat(sprintf("  Analyze Compounds: %s\n", analyze_compounds))
  }

  # No need to generate a filename here, it's done inside the loop
  # output_path <- generate_output_filename(config, args$disease_code, args$output_file)

  if (args$dry_run) {
    cat("\n🧪 DRY RUN MODE - Validation Only\n")
    cat("================================\n")
    cat(sprintf("Output would be saved to dynamically generated files based on template: %s\n", config$output$filename_template))
    cat("All parameters validated successfully ✓\n")
    return(0)
  }

  # The output file is now initialized inside the loop for each compound/category
  # Data is now loaded inside the loop for each measure type

  # Load mapping data separately to determine which categories/compounds to run
  mapping_path <- here(config$data_paths$pesticide_mapping)
  if (!file.exists(mapping_path)) {
    cat(sprintf("❌ Critical error: Pesticide mapping file not found at %s\n", mapping_path))
    return(1)
  }
  mapping_data <- read_csv(mapping_path, show_col_types = FALSE)
  cat(sprintf("✓ Pesticide mapping data loaded successfully (%d records)\n", nrow(mapping_data)))

  if (analyze_compounds) {
    # Extract compound IDs from the argument string (e.g., "compound:1,2,3")
    compound_ids <- gsub("^compound:", "", args$pesticide_category)

    # Get the specific pesticide compounds to analyze based on IDs
    compounds_to_analyze <- get_pesticide_compounds(mapping_data, compound_ids)

    if (nrow(compounds_to_analyze) == 0) {
      cat("❌ No compounds found matching the specified criteria.\n")
      return(1)
    }

    # Calculate the total number of combinations for the entire analysis run
    num_compounds <- nrow(compounds_to_analyze)
    combinations_per_compound <- length(parse_list_argument(args$measure_type)) *
                                 length(parse_list_argument(args$estimate_types)) *
                                 length(parse_list_argument(args$lag_years)) *
                                 length(parse_list_argument(args$model_types))
    total_combinations <- num_compounds * combinations_per_compound

    # Print the initial setup information and table header for the dashboard
    print_analysis_setup(args, config, compounds_to_analyze, total_combinations)
    print_results_table_header()

    # Pre-load all data for all measure types
    cat("\n📂 Loading all data files...\n")
    measure_types <- parse_list_argument(args$measure_type, c('Weight', 'Density'))
    preloaded_data <- list()
    for (measure in measure_types) {
      cat(sprintf("  Loading %s data...\n", measure))
      preloaded_data[[measure]] <- load_all_data(config, args$disease_code, measure)
      if (is.null(preloaded_data[[measure]])) {
        cat(sprintf("  ⚠️  Warning: Failed to load %s data\n", measure))
      }
    }
    cat("✓ Data loading complete\n\n")

    # Initialize the progress bar for the entire run
    pb <- progress_bar$new(
      format = "[:bar] :percent | ETA: :eta | :what",
      total = total_combinations,
      width = 60
    )

    all_failed_combinations <- list()
    total_successful <- 0

    # Loop through each compound selected for analysis
    for (i in 1:num_compounds) {
      compound_info <- compounds_to_analyze[i, ]

      # Announce the start of processing for the current compound
      print_compound_header(i, num_compounds, compound_info$compound_name)

      # Generate a unique, timestamped output filename for this specific compound
      output_path <- generate_output_filename(config, args$disease_code, compound_info$compound_name)
      initialize_output_file(output_path, config)

      # Call the core analysis function for the single compound
      failed_runs <- analyze_single_compound(compound_info, args, config, output_path, pb, preloaded_data)

      # Aggregate results and failures
      all_failed_combinations <- c(all_failed_combinations, failed_runs)
      successful_runs <- combinations_per_compound - length(failed_runs)
      total_successful <- total_successful + successful_runs

      # Print a summary for the completed compound
      print_compound_summary(output_path, successful_runs, combinations_per_compound)
    }

    # After all compounds are processed, print the final summary of the entire run
    print_final_summary(all_failed_combinations, total_successful, total_combinations, analysis_start_time)

    # Return an exit code: 0 for success (at least one combination worked), 1 for total failure
    return(ifelse(total_successful > 0, 0, 1))

  } else {
    # Get pesticide categories to analyze
    categories_to_analyze <- get_pesticide_categories(
      mapping_data, args$pesticide_category)

    if (nrow(categories_to_analyze) == 0) {
      cat("❌ No categories found matching criteria\n")
      return(1)
    }

    cat(sprintf("\n📋 Categories to analyze: %d\n", nrow(categories_to_analyze)))

    # Initialize counters for category analysis
    total_categories <- nrow(categories_to_analyze)
    total_successful_categories <- 0
    all_category_results <- list()

    # Pre-load all data for all measure types
    cat("\n📂 Loading all data files...\n")
    measure_types <- parse_list_argument(args$measure_type, c('Weight', 'Density'))
    preloaded_data <- list()
    for (measure in measure_types) {
      cat(sprintf("  Loading %s data...\n", measure))
      preloaded_data[[measure]] <- load_all_data(config, args$disease_code, measure)
      if (is.null(preloaded_data[[measure]])) {
        cat(sprintf("  ⚠️  Warning: Failed to load %s data\n", measure))
      }
    }
    cat("✓ Data loading complete\n\n")

    # Process each category
    for (i in 1:total_categories) {
      category_info <- categories_to_analyze[i, ]

      cat(sprintf("\n🔄 Processing category %d/%d: %s (ID: %d)\n",
                  i, total_categories, category_info$category, category_info$cat_id))

      # Generate output filename for this category
      output_path <- generate_output_filename(config, args$disease_code, category_info$category)

      # Analyze this category
      category_result <- analyze_single_category(category_info, args, config, output_path, preloaded_data)

      # Store results
      all_category_results[[category_info$category]] <- category_result
      if (category_result$successful > 0) {
        total_successful_categories <- total_successful_categories + 1
      }

      cat(sprintf("  📈 Category %s: %d/%d successful (%.1f%%)\n",
                  category_info$category,
                  category_result$successful,
                  category_result$total,
                  category_result$success_rate))
    }

    # Print final summary for category analysis
    cat("\n📊 FINAL CATEGORY ANALYSIS SUMMARY\n")
    cat(sprintf("  Total categories processed: %d\n", total_categories))
    cat(sprintf("  Successful categories: %d\n", total_successful_categories))
    cat(sprintf("  Success rate: %.1f%%\n",
                ifelse(total_categories > 0, (total_successful_categories/total_categories)*100, 0)))

    # Return success if at least one category was processed successfully
    return(ifelse(total_successful_categories > 0, 0, 1))
  }
}

# Execute main function
if (!interactive()) {
  quit(status = main())
} else {
  main()
}
