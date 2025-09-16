#!/bin/bash
#SBATCH --job-name=WDP_24D_Test
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/public/home/acf4pijnzl/WDP_Results/slurm-%j.out
#SBATCH --error=/public/home/acf4pijnzl/WDP_Results/slurm-%j.err

# WDP 2,4-D Test Analysis Script for SLURM
# This script runs a test analysis on the HPC cluster using SLURM

echo "============================================================"
echo "WDP 2,4-D Test Analysis on SLURM Cluster"
echo "============================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Start time: $(date)"
echo ""

# Load required modules
echo "Loading modules..."
module load R/4.3.0 2>/dev/null || module load R 2>/dev/null || echo "Using system R"

# Check R installation
echo "R version:"
R --version | head -n 1

# Install required packages if not available
echo ""
echo "📦 Checking R packages..."
R --quiet --slave << 'RSCRIPT'
# Function to check and install packages
check_packages <- function() {
    required_packages <- c("INLA", "dplyr", "readr", "yaml", "argparse", "progress")

    # Check which packages are missing
    missing <- required_packages[!required_packages %in% installed.packages()[,"Package"]]

    if (length(missing) > 0) {
        cat("Installing missing packages:", paste(missing, collapse=", "), "\n")

        # Install INLA if missing
        if ("INLA" %in% missing) {
            install.packages("INLA",
                repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"))
        }

        # Install other packages
        other_missing <- missing[missing != "INLA"]
        if (length(other_missing) > 0) {
            install.packages(other_missing, repos = "https://cran.rstudio.com/")
        }
    }

    # Verify all packages can be loaded
    for (pkg in required_packages) {
        if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
            stop(paste("Failed to load package:", pkg))
        }
    }

    cat("✓ All required packages are available\n")
}

check_packages()
RSCRIPT

# Set working directory
cd /public/home/acf4pijnzl/WDP_Analysis

# Create results directory with timestamp
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
RESULTS_DIR="/public/home/acf4pijnzl/WDP_Results/Test_Analysis_${TIMESTAMP}"
mkdir -p "$RESULTS_DIR"

echo "📋 Analysis Parameters:"
echo "----------------------"
echo "  Working directory: $(pwd)"
echo "  Results directory: ${RESULTS_DIR}"
echo "  Output files will be saved to: ${RESULTS_DIR}"
echo ""

# Run test analysis with a single model to verify everything works
echo "🚀 Running test analysis..."
echo "   This will run a single model (M0) with avg estimate for Weight measure"
echo ""

Rscript Code/INLA/BYM_INLA_Production.R \
    --pesticide-category TEST \
    --measure-type Weight \
    --estimate-types avg \
    --lag-years 5 \
    --model-types M0 \
    --disease-code C81-C96 \
    --output-file "${RESULTS_DIR}/Results_Test_M0_Weight_avg_5yr.csv" \
    --verbose 2>&1 | tee "${RESULTS_DIR}/test_analysis.log"

# Check if output was created
if [ -f "${RESULTS_DIR}/Results_Test_M0_Weight_avg_5yr.csv" ]; then
    num_results=$(wc -l < "${RESULTS_DIR}/Results_Test_M0_Weight_avg_5yr.csv")
    echo "   ✓ Test results saved: $((num_results - 1)) records"
else
    echo "   ❌ Test analysis failed - no output file created"
fi

# Generate summary report
echo ""
echo "📊 Generating Summary Report"
echo "============================"

SUMMARY_FILE="${RESULTS_DIR}/Test_Summary_Report.txt"

cat > "$SUMMARY_FILE" << EOF
WDP 2,4-D Test Analysis Summary
==============================
Job ID: $SLURM_JOB_ID
Node: $SLURMD_NODENAME
Generated: $(date)

Analysis Parameters:
- Pesticide Category: TEST
- Measure: Weight
- Estimate: avg
- Lag: 5 years
- Model: M0
- Disease: C81-C96

Files Generated:
EOF

# Add file information
for file in ${RESULTS_DIR}/*; do
    if [ -f "$file" ]; then
        filename=$(basename "$file")
        filesize=$(du -h "$file" | cut -f1)
        echo "- $filename: $filesize" >> "$SUMMARY_FILE"
    fi
done

# Display summary
echo ""
cat "$SUMMARY_FILE"

# Final summary
echo ""
echo "============================================================"
echo "Test Analysis Complete!"
echo "============================================================"
echo "Results directory: ${RESULTS_DIR}"
echo ""
echo "Files created:"
ls -la "${RESULTS_DIR}/"
echo ""
echo "Completed at: $(date)"
echo "============================================================"