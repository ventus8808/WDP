# WDP INLA Analysis - Troubleshooting Guide

## Common Issues and Solutions

### 1. Path-Related Errors

#### Problem: "Configuration file not found" or "Required file not found"

**Symptoms:**
```
Error: Configuration file not found: INLA_Config/analysis_config.yaml
Error: Required file not found: Data/Processed/...
```

**Solution:**
1. Ensure you're running the script from the correct directory:
   ```bash
   cd /path/to/WDP/Code/INLA
   pwd  # Should show: .../WDP/Code/INLA
   ```

2. Verify project structure is intact:
   ```bash
   ls -la INLA_Config/analysis_config.yaml
   ls -la ../../Data/Processed/
   ```

3. Check `here` package installation:
   ```r
   Rscript -e "library(here); print(here())"
   ```

#### Problem: "cannot open the connection" when creating spatial graphs

**Symptoms:**
```
Error in inla.write.graph(): cannot open the connection
```

**Solution:**
1. Ensure temp directories exist and are writable:
   ```bash
   mkdir -p Code/INLA/INLA_Temp/graphs
   chmod 755 Code/INLA/INLA_Temp/graphs
   ```

2. Check disk space:
   ```bash
   df -h .
   ```

3. Verify graph directory configuration in `analysis_config.yaml`:
   ```yaml
   model_fitting:
     spatial:
       graph_dir: "Code/INLA/INLA_Temp/graphs"
   ```

### 2. Data Processing Errors

#### Problem: "rolling_mean" function errors

**Symptoms:**
```
Error in argument: `pesticide_lagged = rolling_mean(.data[["cat33_avg"]], lag_years, align = "right")`
```

**Causes and Solutions:**

1. **Missing data column:**
   - Check if the pesticide column exists in data
   - Verify compound ID is valid in mapping file

2. **Insufficient data for lag calculation:**
   - Ensure data spans more years than lag period
   - Check for data gaps in time series

3. **Data type issues:**
   - Verify numeric data types:
   ```r
   str(pesticide_data$cat33_avg)  # Should be numeric
   sum(is.na(pesticide_data$cat33_avg))  # Check for NAs
   ```

#### Problem: "Insufficient data after processing"

**Symptoms:**
```
Warning: Low record count after merge: 45 records
Error: Insufficient data
```

**Solution:**
1. Check minimum thresholds in config:
   ```yaml
   data_processing:
     min_thresholds:
       records_per_analysis: 100
       counties_per_analysis: 50
   ```

2. Verify data coverage:
   ```r
   # Count available counties
   length(unique(cdc_data$COUNTY_FIPS))
   # Check year range
   range(cdc_data$Year)
   ```

3. Consider relaxing thresholds for exploratory analysis (with caution)

### 3. Model Fitting Issues

#### Problem: INLA model convergence failures

**Symptoms:**
```
Warning: INLA model failed to converge
Error: Model fitting failed
```

**Solutions:**

1. **Check for extreme values:**
   ```r
   summary(model_data$pesticide_lagged)
   # Look for outliers or infinite values
   ```

2. **Verify spatial structure:**
   - Ensure counties are properly connected
   - Check for isolated regions

3. **Simplify model complexity:**
   - Start with M0 (base model)
   - Add covariates gradually

4. **Increase INLA iterations:**
   ```yaml
   model_fitting:
     inla:
       control_compute:
         strategy: "adaptive"
       control_inla:
         strategy: "gaussian"
   ```

### 4. Memory and Performance Issues

#### Problem: "Cannot allocate memory" or slow performance

**Solutions:**

1. **Reduce data size:**
   - Subset to specific years or regions
   - Use fewer model combinations

2. **Optimize INLA settings:**
   ```yaml
   model_fitting:
     inla:
       num_threads: "2:1"  # Reduce threads
   ```

3. **Monitor memory usage:**
   ```bash
   # Check available memory
   free -h
   # Monitor during execution
   top -p $(pgrep -f INLA_Main.R)
   ```

4. **Use SLURM resources effectively:**
   ```bash
   #SBATCH --mem-per-cpu=4G  # Increase memory
   #SBATCH --time=4:00:00    # Increase time limit
   ```

### 5. Package and Environment Issues

#### Problem: Package loading errors

**Symptoms:**
```
Error: there is no package called 'here'
Error in library(INLA): package 'INLA' not found
```

**Solutions:**

1. **Run dependency installer:**
   ```bash
   Rscript INLA_Dependencies/install_packages.R
   ```

2. **For INLA-specific issues, see the complete installation guide:**
   ```
   INLA_Docs/INLA_Installation_Guide.md
   ```

3. **Manual INLA installation:**
   ```r
   install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
   ```

3. **Check R version compatibility:**
   ```r
   R.version.string  # Should be R 4.0+
   ```

4. **Clear package cache if needed:**
   ```r
   remove.packages("INLA")
   # Reinstall
   ```

### 6. SLURM and Server Issues

#### Problem: Job submission failures

**Solutions:**

1. **Check SLURM configuration:**
   ```bash
   sinfo  # Check available partitions
   squeue -u $USER  # Check job queue
   ```

2. **Verify module loading:**
   ```bash
   module avail  # Check available modules
   module list   # Check loaded modules
   ```

3. **Test conda environment:**
   ```bash
   conda activate INLA
   which Rscript
   R --version
   ```

4. **Check file permissions:**
   ```bash
   ls -la INLA_Scripts/run_single_compound.sh
   # Should show execute permissions (-rwxr-xr-x)
   ```

### 7. Output and Results Issues

#### Problem: Empty or corrupted output files

**Solutions:**

1. **Check output directory permissions:**
   ```bash
   ls -ld Result/INLA_Analysis/
   mkdir -p Result/INLA_Analysis
   ```

2. **Verify CSV file format:**
   ```bash
   head -5 Result/INLA_Analysis/Results_*.csv
   ```

3. **Check for partial runs:**
   - Look for incomplete timestamp patterns
   - Verify all model combinations completed

### 8. Debugging Strategies

#### Enable Maximum Logging
```bash
Rscript INLA_Main.R \
  --verbose \
  --disease-code "C81-C96" \
  --pesticide-category "compound:1"
```

#### Use Dry Run Mode
```bash
Rscript INLA_Main.R \
  --dry-run \
  --disease-code "C81-C96" \
  --pesticide-category "compound:1"
```

#### Interactive R Session
```r
# Load the main script environment
source("INLA_Main.R")

# Test individual functions
config <- load_config("INLA_Config/analysis_config.yaml")
data_list <- load_all_data(config, "C81-C96", "Weight")

# Examine data structure
str(data_list)
summary(data_list$pesticide)
```

#### Check System Resources
```bash
# Memory usage
free -h

# Disk space
df -h /tmp
df -h .

# CPU info
lscpu

# Currently running processes
ps aux | grep R
```

### 9. Performance Optimization

#### For Large Datasets
1. **Use data subsetting:**
   ```yaml
   data_processing:
     year_range:
       start: 2010  # Limit years
       end: 2019
   ```

2. **Parallel processing:**
   ```yaml
   model_fitting:
     inla:
       num_threads: "4:1"
   ```

3. **Reduce model complexity:**
   - Start with fewer covariates
   - Use fewer exposure estimates
   - Limit lag periods

#### For Server Environments
1. **Optimize SLURM settings:**
   ```bash
   #SBATCH --cpus-per-task=8
   #SBATCH --mem-per-cpu=3G
   #SBATCH --time=6:00:00
   ```

2. **Use local temp directories:**
   ```bash
   export TMPDIR=/local/scratch/$USER
   mkdir -p $TMPDIR
   ```

### 10. Getting Help

#### Log File Analysis
Always check both stdout and stderr logs:
```bash
# Standard output
cat WDP_INLA_Analysis-123456.log

# Error output  
cat WDP_INLA_Analysis-123456.err
```

#### Contact Information
For persistent issues:
1. Save complete error logs
2. Document the exact command used
3. Note the system environment (R version, OS, etc.)
4. Include data characteristics (size, time range, etc.)

#### Useful Commands for Diagnosis
```bash
# System information
uname -a
R --version
conda --version

# Package information
Rscript -e "sessionInfo()"

# File structure check
find Code/INLA -name "*.R" -exec wc -l {} +
```

Remember: Most issues are related to paths, data availability, or resource constraints. The new structure with `here` package and robust error handling should prevent most common problems.