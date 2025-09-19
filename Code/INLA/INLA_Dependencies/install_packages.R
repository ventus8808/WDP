#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# WDP INLA Environment Setup Script
# Automatically install and manage R package dependencies
# Author: WDP Analysis Team
# Date: 2024

cat("🚀 WDP INLA Environment Setup\n")
cat("============================\n\n")

# Configure CRAN mirror to avoid installation errors
options(repos = c(CRAN = "https://cloud.r-project.org/"))
cat("🌐 CRAN mirror configured: https://cloud.r-project.org/\n\n")

# Required packages for WDP INLA analysis
required_packages <- c(
  "here",           # Robust path management
  "dplyr",          # Data manipulation
  "readr",          # CSV reading
  "yaml",           # Configuration files
  "argparse",       # Command line parsing
  "progress",       # Progress bars
  "INLA"            # Bayesian spatial modeling
)

# Optional packages for enhanced functionality
optional_packages <- c(
  "renv",           # Environment management
  "devtools",       # Development tools
  "testthat",       # Testing framework
  "knitr",          # Documentation
  "rmarkdown"       # Report generation
)

#' Install package if not already installed
#' @param package_name Name of the package to install
#' @param from_cran Whether to install from CRAN (default) or special repo
install_if_missing <- function(package_name, from_cran = TRUE) {
  if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("📦 Installing %s...\n", package_name))
    
    tryCatch({
      if (package_name == "INLA") {
        # INLA requires special installation from dedicated repository
        cat("  🔧 Installing INLA from official repository...\n")
        install.packages("INLA", 
                        repos = c(getOption("repos"), 
                                 INLA = "https://inla.r-inla-download.org/R/stable"), 
                        dep = TRUE)
      } else if (from_cran) {
        install.packages(package_name, dependencies = TRUE)
      }
      
      # Verify installation by loading the package
      if (require(package_name, character.only = TRUE, quietly = TRUE)) {
        cat(sprintf("  ✅ %s installed successfully\n", package_name))
        return(TRUE)
      } else {
        cat(sprintf("  ❌ %s installation verification failed\n", package_name))
        return(FALSE)
      }
    }, error = function(e) {
      cat(sprintf("  ❌ Failed to install %s: %s\n", package_name, e$message))
      return(FALSE)
    })
  } else {
    cat(sprintf("  ✅ %s already installed\n", package_name))
    return(TRUE)
  }
}

#' Check and display package versions
#' @param packages Vector of package names to check
check_package_versions <- function(packages) {
  cat("\n📋 Package Version Report:\n")
  cat("-------------------------\n")
  
  for (pkg in packages) {
    if (require(pkg, character.only = TRUE, quietly = TRUE)) {
      version <- packageVersion(pkg)
      cat(sprintf("  %s: %s\n", pkg, version))
    } else {
      cat(sprintf("  %s: NOT INSTALLED\n", pkg))
    }
  }
}

#' Setup renv for reproducible environments
setup_renv <- function() {
  cat("\n🔒 Setting up renv for reproducible environments...\n")
  
  if (!require("renv", quietly = TRUE)) {
    install.packages("renv")
  }
  
  tryCatch({
    # Initialize renv if not already initialized
    if (!file.exists("renv.lock")) {
      renv::init()
      cat("  ✅ renv initialized\n")
    } else {
      cat("  ✅ renv already initialized\n")
    }
    
    # Take a snapshot of current packages
    renv::snapshot()
    cat("  ✅ Package snapshot created\n")
    
  }, error = function(e) {
    cat(sprintf("  ⚠️  renv setup failed: %s\n", e$message))
  })
}

#' Create a dependency check script
create_dependency_checker <- function() {
  checker_script <- '#!/usr/bin/env Rscript
# Quick dependency check for WDP INLA
required <- c("here", "dplyr", "readr", "yaml", "argparse", "progress", "INLA")
missing <- character(0)

for (pkg in required) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    missing <- c(missing, pkg)
  }
}

if (length(missing) > 0) {
  cat("❌ Missing packages:", paste(missing, collapse = ", "), "\\n")
  cat("Run: Rscript INLA_Dependencies/install_packages.R\\n")
  quit(status = 1)
} else {
  cat("✅ All dependencies satisfied\\n")
  quit(status = 0)
}
'
  
  writeLines(checker_script, "INLA_Dependencies/check_dependencies.R")
  Sys.chmod("INLA_Dependencies/check_dependencies.R", mode = "0755")
  cat("📝 Created dependency checker script\n")
}

# Main installation process
main <- function() {
  cat("🔍 Checking required packages...\n")
  
  success_count <- 0
  total_required <- length(required_packages)
  
  # Install required packages
  for (pkg in required_packages) {
    if (install_if_missing(pkg)) {
      success_count <- success_count + 1
    }
  }
  
  cat(sprintf("\n📊 Installation Summary: %d/%d required packages installed\n", 
              success_count, total_required))
  
  if (success_count == total_required) {
    cat("🎉 All required packages installed successfully!\n")
    
    # Check package versions
    check_package_versions(required_packages)
    
    # Ask about optional packages
    cat("\n💡 Optional packages can enhance functionality. Install them? (y/n): ")
    if (interactive()) {
      response <- readline()
      if (tolower(response) %in% c("y", "yes")) {
        cat("\n📦 Installing optional packages...\n")
        for (pkg in optional_packages) {
          install_if_missing(pkg)
        }
      }
    } else {
      cat("Running in non-interactive mode, skipping optional packages\n")
    }
    
    # Setup renv
    cat("\n🔄 Do you want to setup renv for reproducible environments? (y/n): ")
    if (interactive()) {
      response <- readline()
      if (tolower(response) %in% c("y", "yes")) {
        setup_renv()
      }
    } else {
      cat("Running in non-interactive mode, skipping renv setup\n")
    }
    
    # Create dependency checker
    create_dependency_checker()
    
    cat("\n✨ Environment setup complete!\n")
    cat("📄 Quick start: Rscript INLA_Main.R --help\n")
    cat("🔍 Check deps: Rscript INLA_Dependencies/check_dependencies.R\n")
    
  } else {
    cat("❌ Some required packages failed to install. Please check the errors above.\n")
    return(1)
  }
  
  return(0)
}

# Run main function
if (!interactive()) {
  quit(status = main())
} else {
  main()
}