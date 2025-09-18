# INLA Package Installation Guide for HPC Servers

## Overview

Installing R-INLA (Integrated Nested Laplace Approximation) on HPC servers can be challenging due to complex dependencies and version compatibility issues. This guide provides a proven workflow that combines Conda environment management with R's native package installer to achieve a stable INLA installation.

## Core Strategy

The key insight is **division of labor**:
- **Conda**: Handles system-level dependencies and provides a stable R environment
- **R Package Manager**: Handles specialized packages like INLA from their dedicated repositories

This approach avoids conflicts between different package managers while ensuring version compatibility.

---

## Complete Installation Workflow

### Step 1: Create a Clean Conda Environment with System Dependencies

Create a fresh Conda environment that includes R and system-level packages that are difficult to compile within R:

```bash
# Create new environment with essential system packages
conda create -n INLA \
  r-base \
  r-sf \
  r-terra \
  r-matrix \
  r-dplyr \
  r-readr \
  r-yaml \
  -c conda-forge

# Activate the environment
conda activate INLA
```

**Key Points:**
- We deliberately **exclude** `r-inla` and `r-fmesher` from Conda installation
- Include packages that require complex system libraries (like `r-sf` for geospatial work)
- This provides a stable foundation without version conflicts

### Step 2: Configure R Package Repositories (Critical Step)

This is the most important step. Configure R to use the correct package sources in the right priority order:

```r
# Start R within the Conda environment
R

# Configure package repositories with proper priority
options(repos = c(
  "inlabru-org" = "https://inlabru-org.r-universe.dev",     # Priority 1: Latest fmesher
  "INLA" = "https://inla.r-inla-download.org/R/stable",     # Priority 2: Official INLA
  "CRAN" = "https://cloud.r-project.org"                    # Priority 3: Standard packages
))

# Verify the configuration
getOption("repos")
```

**Why This Works:**
1. **inlabru-org r-universe**: Provides the latest `fmesher` that's compatible with testing INLA
2. **INLA official repository**: Provides the authoritative INLA package
3. **CRAN**: Provides all other standard R packages

This configuration creates a "shopping map" for R to find compatible versions.

### Step 3: Install INLA Using Standard R Commands

With the repositories configured, use R's native installer:

```r
# Install INLA and all dependencies
install.packages("INLA", dependencies = TRUE)

# Verify installation
library(INLA)
inla.version()
```

**What Happens:**
- R automatically downloads the latest compatible `fmesher` from inlabru-org
- R then installs `INLA` from the official repository
- All versions are automatically matched and compatible

### Step 4: Verification and Testing

Test the installation to ensure everything works:

```r
# Load essential packages
library(INLA)
library(sf)
library(dplyr)

# Quick functionality test
test_data <- data.frame(
  y = rpois(100, 10),
  x = rnorm(100),
  idx = 1:100
)

# Simple INLA model test
formula <- y ~ x + f(idx, model="rw1")
result <- inla(formula, data=test_data, family="poisson")

# Check if model ran successfully
summary(result)
```

---

## Common Issues and Solutions

### Problem: Version Mismatch Between INLA and fmesher

**Symptoms:**
```
Error in fm_identical_CRS(...)
Error: function 'fm_...' not found
```

**Solution:**
- Ensure you're using the repository configuration from Step 2
- The inlabru-org r-universe typically has the most up-to-date fmesher

### Problem: Compilation Errors

**Symptoms:**
```
ERROR: compilation failed for package 'INLA'
```

**Solution:**
- Use Conda to install system dependencies first (Step 1)
- Ensure you have proper compiler tools loaded:
```bash
module load compiler/dtk/23.10  # or appropriate for your system
```

### Problem: Network/Download Issues

**Symptoms:**
```
Warning: unable to access index for repository
```

**Solution:**
- Verify internet connectivity from compute nodes
- Try alternative INLA repository URLs:
```r
options(repos = c(
  "INLA" = "https://inla.r-inla-download.org/R/testing",  # Use testing if stable fails
  "CRAN" = "https://cloud.r-project.org"
))
```

---

## Environment Management Best Practices

### Save Your Configuration

Create a startup script to preserve your configuration:

```r
# Create .Rprofile in your project directory
cat('
# INLA-optimized repository configuration
local({
  repos <- c(
    "inlabru-org" = "https://inlabru-org.r-universe.dev",
    "INLA" = "https://inla.r-inla-download.org/R/stable", 
    "CRAN" = "https://cloud.r-project.org"
  )
  options(repos = repos)
})
', file = "~/.Rprofile")
```

### Document Your Environment

Record successful configurations:

```bash
# Save package versions
Rscript -e "sessionInfo()" > environment_info.txt

# Save Conda environment
conda env export > environment.yml
```

### Automated Installation Script

Create a reusable installation script:

```bash
#!/bin/bash
# inla_install.sh

echo "Setting up INLA environment..."

# Create Conda environment
conda create -n INLA r-base r-sf r-terra r-matrix -c conda-forge -y
conda activate INLA

# Configure and install INLA
Rscript -e '
options(repos = c(
  "inlabru-org" = "https://inlabru-org.r-universe.dev",
  "INLA" = "https://inla.r-inla-download.org/R/stable",
  "CRAN" = "https://cloud.r-project.org"
))
install.packages("INLA", dependencies = TRUE)
library(INLA)
cat("INLA version:", inla.version(), "\n")
'

echo "INLA installation complete!"
```

---

## Integration with WDP INLA System

### Using with the Restructured WDP System

The WDP INLA analysis system now includes automated dependency management:

```bash
# Use the automated installer
Rscript INLA_Dependencies/install_packages.R

# Or check existing installation
Rscript INLA_Dependencies/check_dependencies.R
```

### SLURM Integration

The optimized SLURM script handles environment activation:

```bash
# Submit job with proper environment
sbatch INLA_Scripts/run_single_compound.sh C81-C96 1
```

---

## Troubleshooting Checklist

When INLA installation fails, check these in order:

1. **Environment**: Is Conda activated correctly?
2. **Repositories**: Are the package sources configured properly?
3. **Network**: Can you reach the INLA repositories?
4. **Dependencies**: Are system libraries available through Conda?
5. **Permissions**: Can you write to the package installation directory?
6. **Versions**: Are you using compatible R and package versions?

### Getting Help

For persistent issues:
1. Save complete error logs
2. Document your system configuration (OS, R version, etc.)
3. Note the exact commands used
4. Check the INLA user forum: https://groups.google.com/g/r-inla-discussion-group

---

## Summary

**The Golden Rule**: Let Conda handle the foundation, let R handle the specialization.

This approach has proven successful across multiple HPC environments and resolves the most common INLA installation challenges. The key is understanding that INLA and its ecosystem move faster than traditional package managers, so using the specialized repositories is essential for success.