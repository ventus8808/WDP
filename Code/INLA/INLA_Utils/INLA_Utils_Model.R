#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# WDP BYM INLA Model Fitting Utilities
# Functions for spatial structure creation, model fitting, and validation
# Author: WDP Analysis Team
# Date: 2024

# Load required libraries
suppressMessages({
  if (!require(here, quietly = TRUE)) {
    install.packages("here")
    library(here, quietly = TRUE)
  }
})

#' Create spatial adjacency structure for INLA analysis
#' @param adjacency_data County adjacency data frame
#' @param counties_in_data Vector of county FIPS codes in analysis
#' @param category_id Unique identifier for graph file naming
#' @param config Configuration list
#' @return INLA graph object and file path
create_spatial_structure <- function(adjacency_data, counties_in_data, category_id, config) {
  cat("  🗺️  Creating spatial structure...\n")

  suppressWarnings({
    # <<<<<<<  新增诊断代码开始  >>>>>>>
    cat("\n--- [DEBUG] Checking input parameters in create_spatial_structure ---\n")
    cat(sprintf("Length of counties_in_data: %d\n", length(counties_in_data)))
    cat(sprintf("Class of counties_in_data: %s\n", class(counties_in_data)))
    cat("First 10 county FIPS in counties_in_data:\n")
    if(length(counties_in_data) > 0) {
      print(head(counties_in_data, 10))
    }
    cat(sprintf("NROW of adjacency_data: %d\n", nrow(adjacency_data)))
    cat(sprintf("Class of adjacency_data: %s\n", class(adjacency_data)))
    if(nrow(adjacency_data) > 0) {
      cat("First few rows of adjacency_data:\n")
      print(head(adjacency_data, 3))
    }
    cat(sprintf("Category ID: %s\n", category_id))
    cat("--- [DEBUG] End of input parameters check ---\n\n")
    # <<<<<<<  新增诊断代码结束  >>>>>>>
    
    # Sort counties for consistent indexing
    cat("  📊 About to sort counties...\n")
    region_ids <- sort(unique(as.character(counties_in_data)))
    cat(sprintf("  📊 Sorted counties: %d unique regions\n", length(region_ids)))
    n_regions <- length(region_ids)
    cat(sprintf("  📊 Number of regions: %d\n", n_regions))

    # Create FIPS to index mapping
    cat("  📊 Creating FIPS to index mapping...\n")
    fips_to_index <- setNames(1:n_regions, region_ids)
    cat("  📊 FIPS to index mapping created\n")

    # Initialize adjacency matrix
    cat("  📊 Initializing adjacency matrix...\n")
    adj_matrix <- matrix(0, nrow = n_regions, ncol = n_regions)
    cat(sprintf("  📊 Adjacency matrix initialized with dimensions %d x %d\n", 
                nrow(adj_matrix), ncol(adj_matrix)))

    # Filter adjacency data to analysis counties only
    cat("  📊 Filtering adjacency data...\n")
    adj_filtered <- adjacency_data %>%
      filter(
        county_from %in% region_ids,
        county_to %in% region_ids
      )
    cat(sprintf("  📊 Filtered adjacency data has %d rows\n", nrow(adj_filtered)))

    # Fill adjacency matrix
    if (nrow(adj_filtered) > 0) {
      from_indices <- fips_to_index[as.character(adj_filtered$county_from)]
      to_indices <- fips_to_index[as.character(adj_filtered$county_to)]

      # Create adjacency pairs
      adj_pairs <- cbind(from_indices, to_indices)

      # Set adjacency (symmetric)
      adj_matrix[adj_pairs] <- 1
      adj_matrix[adj_pairs[, c(2, 1)]] <- 1
    }
    
    # <<<<<<<  新增诊断代码开始  >>>>>>>
    cat(sprintf("  📊 Adjacency matrix created with dimensions %d x %d\n", 
                nrow(adj_matrix), ncol(adj_matrix)))
    cat(sprintf("  📊 Number of connections in adjacency matrix: %d\n", sum(adj_matrix)))
    # <<<<<<<  新增诊断代码结束  >>>>>>>

    # Create INLA graph
    # <<<<<<<  新增诊断代码开始  >>>>>>>
    cat(sprintf("  📊 About to create INLA graph with adj_matrix of size %d x %d\n", 
                nrow(adj_matrix), ncol(adj_matrix)))
    # 检查邻接矩阵是否为空
    if(nrow(adj_matrix) == 0 || ncol(adj_matrix) == 0) {
      cat("  ❌ Error: Adjacency matrix is empty!\n")
      stop("Adjacency matrix is empty")
    }
    # <<<<<<<  新增诊断代码结束  >>>>>>>
    
    inla_graph <- tryCatch({
      inla.read.graph(adj_matrix)
    }, error = function(e) {
      cat(sprintf("  ❌ Error in inla.read.graph: %s\n", e$message))
      cat("  📊 Adjacency matrix content:\n")
      print(adj_matrix[1:min(5, nrow(adj_matrix)), 1:min(5, ncol(adj_matrix))])
      stop(e$message)
    })
    cat("  📊 INLA graph created successfully\n")

    # Create unique graph file for this analysis
    graph_filename <- sprintf("%s_%s_%d.graph",
                              config$model_fitting$spatial$graph_prefix,
                              "cat", category_id)
    # Resolve graph_dir: if absolute path, use as-is; else resolve relative to project root
    graph_dir <- config$model_fitting$spatial$graph_dir
    if (!dir.exists(graph_dir)) {
      # If it doesn't exist as absolute, try relative to project root
      candidate <- here(graph_dir)
      if (dir.exists(candidate)) {
        graph_dir <- candidate
      }
    }
    if (!dir.exists(graph_dir)) {
      dir.create(graph_dir, recursive = TRUE, showWarnings = FALSE)
    }
    graph_filepath <- file.path(graph_dir, graph_filename)

    # Write graph file
    cat(sprintf("  📊 About to write graph file to: %s\n", graph_filepath))
    # 确保目录存在
    graph_dir <- dirname(graph_filepath)
    if (!dir.exists(graph_dir)) {
      dir.create(graph_dir, recursive = TRUE, showWarnings = FALSE)
      cat(sprintf("  📊 Created directory: %s\n", graph_dir))
    }
    
    inla.write.graph(inla_graph, filename = graph_filepath)
    cat(sprintf("  🗂️  Graph file path: %s (will write/read)\n", graph_filepath))
    cat("  📊 Graph file written successfully\n")

    # Also ensure a copy exists inside INLA working.directory to avoid cross-dir issues
    wd_dir <- tryCatch(inla.getOption("working.directory"), error = function(e) NA)
    if (!is.na(wd_dir) && nzchar(wd_dir)) {
      if (!dir.exists(wd_dir)) dir.create(wd_dir, recursive = TRUE, showWarnings = FALSE)
      wd_graph_path <- file.path(wd_dir, basename(graph_filepath))
      # Copy and ensure permissions
      ok_copy <- tryCatch({
        file.copy(graph_filepath, wd_graph_path, overwrite = TRUE)
      }, error = function(e) FALSE)
      if (ok_copy && file.exists(wd_graph_path)) {
        Sys.chmod(wd_graph_path, mode = "0644", use_umask = TRUE)
        graph_filepath <- wd_graph_path
        cat(sprintf("  📎 Graph also placed in working.dir: %s\n", wd_graph_path))
      } else {
        cat("  ⚠️  Failed to copy graph into working.directory; will use original path.\n")
      }
    }

    # Calculate connectivity statistics
    n_connections <- sum(adj_matrix) / 2  # Divide by 2 for symmetric matrix
    connectivity_rate <- n_connections / (n_regions * (n_regions - 1) / 2)

    cat(sprintf("  ✓ Spatial graph: %d counties, %d connections (%.2f%% connectivity)\n",
                n_regions, n_connections, connectivity_rate * 100))

    # <<<<<<<  新增诊断代码开始  >>>>>>>
    cat("  📊 About to return spatial structure...\n")
    result <- list(
      graph = inla_graph,
      filepath = graph_filepath,
      n_regions = n_regions,
      n_connections = n_connections,
      connectivity_rate = connectivity_rate
    )
    cat("  📊 Spatial structure created successfully\n")
    return(result)
    # <<<<<<<  新增诊断代码结束  >>>>>>>
  })
}

#' Build model formula based on model type and available covariates
#' @param model_type Model identifier (M0, M1, M2, M3)
#' @param model_data Data frame to check for available covariates
#' @param spatial_graph_path Path to spatial graph file
#' @param config Configuration list
#' @return Formula object for INLA
build_model_formula <- function(model_type, model_data, spatial_graph_path, config) {
  cat(sprintf("  🔧 Building %s model formula...\n", model_type))

  # Determine if we should use non-linear model
  use_nonlinear <- if (!is.null(config$model_fitting$nonlinear$enabled)) {
    config$model_fitting$nonlinear$enabled && model_type %in% config$model_fitting$nonlinear$model_types
  } else {
    FALSE
  }

  # Base formula components
  if (use_nonlinear && "pesticide_binned_idx" %in% names(model_data)) {
    # Non-linear dose-response using random walk on binned exposure
    base_components <- c(
      "1",  # Intercept
      sprintf("f(pesticide_binned_idx, model = 'rw2', hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))")
    )
    cat("    Using non-linear (RW2) dose-response model\n")
  } else {
    # Linear dose-response model with log-transformed standardized exposure
    base_components <- c(
      "1",  # Intercept
      "pesticide_log_std"  # Linear dose-response on log scale
    )
    cat("    Using linear dose-response model\n")
  }

  # 极简模型：不加空间和时间项，只保留固定效应和协变量
  # spatial_component <- "f(county_idx, model = 'iid')" # 注释掉
  # temporal_component <- sprintf("f(Year, model = '%s')", config$model_fitting$temporal$model_type) # 注释掉
  # Get model配置
  model_config <- config$analysis$models[[model_type]]
  if (is.null(model_config)) {
    stop(sprintf("Unknown model type: %s", model_type))
  }

  # Check availability of covariates and add to formula
  covariate_components <- c()
  available_covariates <- names(model_data)

  for (covariate in model_config$covariates) {
    # Map to standardized column names
    std_col_name <- switch(covariate,
      "SVI_PCA" = "SVI_std",
      "Climate_Factor_1" = "Climate1_std",
      "Climate_Factor_2" = "Climate2_std",
      covariate  # Default to original name
    )

    if (std_col_name %in% available_covariates) {
      # Check if covariate has variation (not all NA)
      if (!all(is.na(model_data[[std_col_name]]))) {
        covariate_components <- c(covariate_components, std_col_name)
      } else {
        warning(sprintf("Covariate %s is all NA, excluding from model", std_col_name))
      }
    } else {
      warning(sprintf("Covariate %s not found in data, excluding from model", std_col_name))
    }
  }

  # Check if required exposure variables exist in data
  if (use_nonlinear) {
    if (!"pesticide_binned_idx" %in% names(model_data)) {
      stop("Non-linear exposure variable 'pesticide_binned_idx' not found in model data. Ensure prepare_model_data creates this variable.")
    }
  } else {
    if (!"pesticide_log_std" %in% names(model_data)) {
      stop("Continuous exposure variable 'pesticide_log_std' not found in model data. Ensure prepare_model_data creates this variable.")
    }
  }

  # Combine all components
  formula_components <- c(
    base_components,
    covariate_components
    # spatial_component, # 注释掉
    # temporal_component # 注释掉
  )

  # Build formula string
  formula_string <- paste("Deaths ~", paste(formula_components, collapse = " + "))

  # Convert to formula object
  formula_obj <- as.formula(formula_string)
  # Attach graph path for diagnostics and use basename to avoid cross-dir issues
  attr(formula_obj, "_graph_path") <- spatial_graph_path
  attr(formula_obj, "_graph_basename") <- basename(spatial_graph_path)

  cat(sprintf("  ✓ Formula: %s\n", deparse(formula_obj)[1]))
  cat(sprintf("    Covariates included: %s\n",
              ifelse(length(covariate_components) > 0,
                     paste(covariate_components, collapse = ", "),
                     "None")))

  return(formula_obj)
}

#' Fit INLA model with comprehensive error handling
#' @param formula Model formula
#' @param model_data Data frame for analysis
#' @param config Configuration list
#' @return INLA model object or NULL if failed
fit_inla_model <- function(formula, model_data, config) {
  # 打印 graph 文件内容和 md5，便于诊断格式问题
  graph_path <- attr(formula, "_graph_path")
  if (!is.null(graph_path) && file.exists(graph_path)) {
    cat("    [diag] Graph file md5:", tryCatch(tools::md5sum(graph_path), error=function(e) "md5 error"), "\n")
    cat("    [diag] Graph file head:\n")
    cat(paste(tryCatch(readLines(graph_path, n=10), error=function(e) "read error"), collapse="\n"), "\n")
  }
  cat("  🎯 Fitting INLA model...\n")
  # Diagnostics: show and verify working directory and graph file accessibility
  wd_opt <- tryCatch(inla.getOption("working.directory"), error = function(e) NA)
  if (is.na(wd_opt) || !nzchar(wd_opt)) wd_opt <- tempdir()
  if (!dir.exists(wd_opt)) {
    dir.create(wd_opt, recursive = TRUE, showWarnings = FALSE)
  }
  can_write_wd <- tryCatch({
    tf <- file.path(wd_opt, sprintf("_wdp_inla_probe_%s", as.integer(Sys.time())))
    ok <- TRUE
    writeLines("ok", tf)
    ok <- ok && file.exists(tf)
    unlink(tf, force = TRUE)
    ok
  }, error = function(e) FALSE)
  cat(sprintf("    [diag] INLA working.directory=%s, writable=%s\n", wd_opt, as.character(can_write_wd)))

  # Extract graph path from formula attribute (set in build_model_formula)
  graph_path <- attr(formula, "_graph_path")
  if (!is.null(graph_path)) {
    cat(sprintf("    [diag] Graph file exists=%s at %s\n", as.character(file.exists(graph_path)), graph_path))
    # List files in graph directory for visibility
    gdir <- dirname(graph_path)
    if (dir.exists(gdir)) {
          files <- tryCatch(list.files(gdir, pattern = "\\.graph$", full.names = TRUE), error = function(e) character(0))
      cat(sprintf("    [diag] Graph dir: %s, count=%d\n", gdir, length(files)))
    }
  } else {
    cat("    [diag] Graph path attribute not found on formula\n")
  }

  # Show current working directory from R side and ensure it's not on NFS
  cat(sprintf("    [diag] R getwd()=%s\n", getwd()))
  cat(sprintf("    [diag] R tempdir()=%s\n", tempdir()))
  # Print a shallow tree of working.directory for debugging
  wd_ls <- tryCatch(list.files(wd_opt, all.files = TRUE), error = function(e) character(0))
  cat(sprintf("    [diag] workdir file count=%d\n", length(wd_ls)))
  
  # Pre-flight checks
  if (nrow(model_data) < 100) {
    cat("  ⚠️ Warning: Very small dataset (< 100 observations)\n")
  }
  
  # Check for data issues that could cause INLA to fail
  numeric_vars <- c("Deaths", "log_expected", "pesticide_log_std")
  for (var in numeric_vars) {
    if (var %in% names(model_data)) {
      if (any(is.infinite(model_data[[var]]) | is.nan(model_data[[var]]))) {
        cat(sprintf("  ❌ Found infinite/NaN values in %s\n", var))
        return(NULL)
      }
    }
  }

  # Use conservative INLA settings
  control_compute <- list(
    dic = TRUE,
    waic = FALSE,  # Disable WAIC to reduce computation
    cpo = FALSE    # Disable CPO to reduce computation
  )
  
  control_predictor <- list(compute = FALSE)  # Disable predictor computation
  
  # Conservative INLA control
  control_inla <- list(
    strategy = "gaussian",  # Use Gaussian approximation for stability
    int.strategy = "eb"     # Use empirical Bayes for hyperparameters
  )

  # Additional convergence controls if specified
  if (!is.null(config$quality_control$convergence$max_iterations)) {
    control_inla$max.iter <- min(config$quality_control$convergence$max_iterations, 100)  # Cap iterations
  }

  # 暂时注释sink()，让所有INLA/C端错误直接输出到日志

  model <- tryCatch({
    # Ensure clean environment for INLA
    gc()  # Garbage collection before model fitting
    
    # Ensure we run in the INLA working directory so C-side resolves relative paths properly
    owd <- getwd()
    setwd(wd_opt)
    cat(sprintf("    [diag] setwd to working.directory: %s (old=%s)\n", wd_opt, owd))

    # Replace graph path in formula with basename to ensure INLA reads from cwd
    gbase <- attr(formula, "_graph_basename")
    if (!is.null(gbase)) {
      # Rebuild formula string swapping the graph path occurrence, then re-parse
      f_str <- deparse(formula)
      f_str <- gsub("graph = '.*?'", sprintf("graph = '%s'", gbase), f_str)
      formula <- as.formula(f_str)
    }

    # The main INLA call with conservative settings
    result <- inla(
      formula = formula,
      data = model_data,
      family = "poisson",
      offset = model_data$log_expected,
      control.compute = control_compute,
      control.predictor = control_predictor,
      control.inla = control_inla,
      verbose = TRUE,
      keep = TRUE,  # Keep intermediate results for stability/diagnostics
      working.directory = wd_opt  # Use INLA global working directory (node-local)
    )
    
    # Clean up after INLA
    gc()
  setwd(owd)
  cat(sprintf("    [diag] restored setwd to: %s\n", owd))
    result
    
  }, error = function(e) {
    cat(sprintf("  ✗ Model fitting error: %s\n", e$message))
    
    # Additional diagnostic information
    if (grepl("file.exists", e$message)) {
      cat("  💡 This appears to be an INLA temporary file issue\n")
      cat("  💡 Suggestions:\n")
      cat("    - Check disk space and permissions\n")
      cat("    - Try reducing model complexity\n")
      cat("    - Verify INLA installation\n")
    }
    
    return(NULL)
  })

  # 暂不恢复sink

  # Log detailed error information if debugging
  if (is.null(model) && config$logging$level %in% c("DEBUG", "INFO")) {
      cat(sprintf("    Formula: %s\n", deparse(formula)[1]))
      cat(sprintf("    Data dimensions: %d x %d\n", nrow(model_data), ncol(model_data)))
      cat(sprintf("    Missing values check:\n"))

      key_vars <- c("Deaths", "log_expected", "exposure_group", "county_idx", "Year")
      for (var in key_vars) {
        if (var %in% names(model_data)) {
          n_missing <- sum(is.na(model_data[[var]]))
          cat(sprintf("      %s: %d missing\n", var, n_missing))
        }
      }
  }

  # Validate model fit
  if (!is.null(model)) {
    validation_result <- validate_model_fit(model, config)
    if (!validation_result$valid) {
      cat(sprintf("  ✗ Model validation failed: %s\n", validation_result$message))
      return(NULL)
    }
    cat("  ✓ Model fitted successfully\n")
  }

  return(model)
}

#' Validate INLA model fit quality and convergence
#' @param model INLA model object
#' @param config Configuration list
#' @return List with validation results
validate_model_fit <- function(model, config) {

  if (is.null(model)) {
    return(list(valid = FALSE, message = "Model is NULL"))
  }

  # Check if model has converged
  if (!is.null(model$mode) && any(!model$mode$mode.status %in% c(0, 1))) {
    return(list(valid = FALSE, message = "Model did not converge"))
  }

  # Check for reasonable model diagnostics
  if (config$quality_control$validation$rr_bounds[1] > 0) {
    # Only check if we have fixed effects
    if (!is.null(model$summary.fixed) && nrow(model$summary.fixed) > 0) {
      # Check if any coefficient produces unreasonable RR
      fixed_effects <- model$summary.fixed
      rr_values <- exp(fixed_effects$mean)

      rr_bounds <- config$quality_control$validation$rr_bounds
      if (any(rr_values < rr_bounds[1] | rr_values > rr_bounds[2], na.rm = TRUE)) {
        extreme_vars <- rownames(fixed_effects)[
          rr_values < rr_bounds[1] | rr_values > rr_bounds[2]
        ]
        warning(sprintf("Extreme RR values detected for: %s",
                       paste(extreme_vars, collapse = ", ")))
      }
    }
  }

  # Check model information criteria
  if (!is.null(model$dic) && !is.finite(model$dic$dic)) {
    return(list(valid = FALSE, message = "Invalid DIC value"))
  }

  if (!is.null(model$waic) && !is.finite(model$waic$waic)) {
    return(list(valid = FALSE, message = "Invalid WAIC value"))
  }

  return(list(valid = TRUE, message = "Model validation passed"))
}

#' Extract model diagnostics and information criteria
#' @param model INLA model object
#' @param config Configuration list
#' @return List of diagnostic values
get_model_diagnostics <- function(model, config) {
  if (is.null(model)) {
    return(list(
      dic = NA,
      waic = NA,
      cpo_failure_rate = NA,
      convergence_status = "FAILED",
      n_fixed_effects = 0,
      n_hyperparameters = 0
    ))
  }

  # Extract information criteria
  dic_value <- if (!is.null(model$dic)) model$dic$dic else NA
  waic_value <- if (!is.null(model$waic)) model$waic$waic else NA

  # Calculate CPO failure rate if available
  cpo_failure_rate <- NA
  if (!is.null(model$cpo) && !is.null(model$cpo$failure)) {
    cpo_failure_rate <- mean(model$cpo$failure, na.rm = TRUE)
  }

  # Determine convergence status
  convergence_status <- "SUCCESS"
  if (!is.null(model$mode) && any(!model$mode$mode.status %in% c(0, 1))) {
    convergence_status <- "CONVERGENCE_WARNING"
  }
  if (is.na(dic_value) && is.na(waic_value)) {
    convergence_status <- "FAILED"
  }

  # Count model components
  n_fixed_effects <- if (!is.null(model$summary.fixed)) nrow(model$summary.fixed) else 0
  n_hyperparameters <- if (!is.null(model$summary.hyperpar)) nrow(model$summary.hyperpar) else 0

  return(list(
    dic = dic_value,
    waic = waic_value,
    cpo_failure_rate = cpo_failure_rate,
    convergence_status = convergence_status,
    n_fixed_effects = n_fixed_effects,
    n_hyperparameters = n_hyperparameters
  ))
}

#' Extract spatial and temporal random effects summaries
#' @param model INLA model object
#' @param config Configuration list
#' @return List with random effect summaries
get_random_effects_summary <- function(model, config) {
  if (is.null(model) || is.null(model$summary.random)) {
    return(list(
      spatial_variance = NA,
      spatial_range = NA,
      temporal_variance = NA,
      temporal_precision = NA
    ))
  }

  random_effects <- list(
    spatial_variance = NA,
    spatial_range = NA,
    temporal_variance = NA,
    temporal_precision = NA
  )

  # Extract spatial effects (BYM2 model)
  if ("county_idx" %in% names(model$summary.random)) {
    spatial_summary <- model$summary.random$county_idx
    if (!is.null(spatial_summary)) {
      random_effects$spatial_variance <- var(spatial_summary$mean, na.rm = TRUE)
    }
  }

  # Extract temporal effects (RW1 model)
  if ("Year" %in% names(model$summary.random)) {
    temporal_summary <- model$summary.random$Year
    if (!is.null(temporal_summary)) {
      random_effects$temporal_variance <- var(temporal_summary$mean, na.rm = TRUE)
    }
  }

  # Extract hyperparameters if available
  if (!is.null(model$summary.hyperpar)) {
    hyperpar <- model$summary.hyperpar

    # BYM2 hyperparameters
    if ("Precision for county_idx" %in% rownames(hyperpar)) {
      random_effects$spatial_precision <- hyperpar["Precision for county_idx", "mean"]
    }

    # Temporal hyperparameters
    if ("Precision for Year" %in% rownames(hyperpar)) {
      random_effects$temporal_precision <- hyperpar["Precision for Year", "mean"]
    }
  }

  return(random_effects)
}

#' Clean up temporary spatial graph files
#' @param graph_filepath Path to graph file to remove
cleanup_spatial_files <- function(graph_filepath) {
  if (file.exists(graph_filepath)) {
    file.remove(graph_filepath)
    cat(sprintf("  🧹 Cleaned up spatial graph: %s\n", basename(graph_filepath)))
  }
}

#' Perform model comparison using information criteria
#' @param model_list List of fitted INLA models
#' @param model_names Names of the models
#' @param config Configuration list
#' @return Data frame with model comparison results
compare_models <- function(model_list, model_names, config) {
  if (length(model_list) != length(model_names)) {
    stop("Number of models and names must match")
  }

  comparison_results <- data.frame(
    Model = model_names,
    DIC = NA,
    WAIC = NA,
    Delta_DIC = NA,
    Delta_WAIC = NA,
    Best_DIC = FALSE,
    Best_WAIC = FALSE,
    stringsAsFactors = FALSE
  )

  # Extract DIC and WAIC values
  for (i in seq_along(model_list)) {
    model <- model_list[[i]]
    if (!is.null(model)) {
      comparison_results$DIC[i] <- if (!is.null(model$dic)) model$dic$dic else NA
      comparison_results$WAIC[i] <- if (!is.null(model$waic)) model$waic$waic else NA
    }
  }

  # Calculate deltas (difference from best model)
  if (any(!is.na(comparison_results$DIC))) {
    min_dic <- min(comparison_results$DIC, na.rm = TRUE)
    comparison_results$Delta_DIC <- comparison_results$DIC - min_dic
    comparison_results$Best_DIC <- comparison_results$DIC == min_dic
  }

  if (any(!is.na(comparison_results$WAIC))) {
    min_waic <- min(comparison_results$WAIC, na.rm = TRUE)
    comparison_results$Delta_WAIC <- comparison_results$WAIC - min_waic
    comparison_results$Best_WAIC <- comparison_results$WAIC == min_waic
  }

  return(comparison_results)
}

cat("✓ Model fitting utilities loaded successfully\n")
