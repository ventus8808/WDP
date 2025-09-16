#!/bin/bash
#SBATCH --job-name=WDP_24D_Full
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=/public/home/acf4pijnzl/WDP_Results/slurm-%j.out
#SBATCH --error=/public/home/acf4pijnzl/WDP_Results/slurm-%j.err

# WDP Comprehensive 2,4-D Analysis Script for SLURM
# This script runs the full 2,4-D analysis on the HPC cluster using SLURM

echo "============================================================"
echo "WDP Comprehensive 2,4-D Analysis on SLURM Cluster"
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

# Set analysis parameters
COMPOUND_ID="2"
COMPOUND_NAME="2,4-D"
DISEASE_CODE="C81-C96"
MEASURE_TYPES="Weight,Density"
ESTIMATE_TYPES="min,avg,max"
LAG_YEARS="5"
MODEL_TYPES="M0,M1,M2,M3"

# Create results directory with timestamp
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
RESULTS_DIR="/public/home/acf4pijnzl/WDP_Results/24D_Full_Analysis_${TIMESTAMP}"
mkdir -p "$RESULTS_DIR"

echo "📋 Analysis Parameters:"
echo "----------------------"
echo "  Compound: ${COMPOUND_NAME} (ID: ${COMPOUND_ID})"
echo "  Disease: ${DISEASE_CODE}"
echo "  Measures: Weight, Density"
echo "  Estimates: min, avg, max"
echo "  Lag: ${LAG_YEARS} years"
echo "  Models: M0, M1, M2, M3"
echo "  Working directory: $(pwd)"
echo "  Results directory: ${RESULTS_DIR}"
echo ""

# Create temporary config directory
CONFIG_DIR="${RESULTS_DIR}/configs"
mkdir -p "$CONFIG_DIR"

# Function to create config files
create_configs() {
    echo "📝 Creating configuration files..."

    # Copy original config
    cp Code/INLA/config/analysis_config.yaml "$CONFIG_DIR/config_base.yaml"

    # Create linear config
    cp "$CONFIG_DIR/config_base.yaml" "$CONFIG_DIR/config_linear.yaml"

    # Create non-linear config for all models
    sed 's/enabled: false # Set to true to enable non-linear models/enabled: true # Set to true to enable non-linear models/' \
        "$CONFIG_DIR/config_base.yaml" > "$CONFIG_DIR/config_nonlinear_all.yaml"
    sed -i 's/model_types: \["M2", "M3"\]/model_types: ["M0", "M1", "M2", "M3"]/' "$CONFIG_DIR/config_nonlinear_all.yaml"

    # Create non-linear config for M2,M3 only
    sed 's/enabled: false # Set to true to enable non-linear models/enabled: true # Set to true to enable non-linear models/' \
        "$CONFIG_DIR/config_base.yaml" > "$CONFIG_DIR/config_nonlinear_m2m3.yaml"

    echo "✓ Configuration files created"
}

# Function to run analysis
run_analysis() {
    local config_name=$1
    local output_suffix=$2
    local description=$3

    echo ""
    echo "🚀 Running: ${description}"
    echo "   Config: ${config_name}"
    echo "   Output: Results_C81-C96_2,4-D_${output_suffix}.csv"

    # Run the analysis
    Rscript Code/INLA/BYM_INLA_Production.R \
        --config "${CONFIG_DIR}/${config_name}" \
        --pesticide-category compound:${COMPOUND_ID} \
        --measure-type ${MEASURE_TYPES} \
        --estimate-types ${ESTIMATE_TYPES} \
        --lag-years ${LAG_YEARS} \
        --model-types ${MODEL_TYPES} \
        --disease-code ${DISEASE_CODE} \
        --output-file ${RESULTS_DIR}/Results_C81-C96_2,4-D_${output_suffix}.csv \
        --verbose 2>&1 | tee "${RESULTS_DIR}/analysis_${output_suffix}.log"

    # Check if output was created
    if [ -f "${RESULTS_DIR}/Results_C81-C96_2,4-D_${output_suffix}.csv" ]; then
        local num_results=$(wc -l < "${RESULTS_DIR}/Results_C81-C96_2,4-D_${output_suffix}.csv")
        echo "   ✓ Results saved: $((num_results - 1)) records"
    else
        echo "   ❌ Analysis failed - no output file created"
    fi
}

# Create configurations
create_configs

# Run analyses
echo ""
echo "🔄 Starting Analysis Runs"
echo "========================"

# 1. Linear models only
run_analysis "config_linear.yaml" \
    "Linear_${TIMESTAMP}" \
    "Linear dose-response models (all combinations)"

# 2. Non-linear models for M2,M3 (mixed approach)
run_analysis "config_nonlinear_m2m3.yaml" \
    "Mixed_${TIMESTAMP}" \
    "Mixed models (Linear: M0,M1 | Non-linear: M2,M3)"

# 3. Non-linear models for all (full non-linear)
run_analysis "config_nonlinear_all.yaml" \
    "Nonlinear_${TIMESTAMP}" \
    "Non-linear dose-response models (all combinations)"

# Generate summary report
echo ""
echo "📊 Generating Summary Report"
echo "============================"

SUMMARY_FILE="${RESULTS_DIR}/Summary_Report.txt"

cat > "$SUMMARY_FILE" << EOF
WDP 2,4-D Comprehensive Analysis Summary
========================================
Job ID: $SLURM_JOB_ID
Node: $SLURMD_NODENAME
Generated: $(date)

Analysis Parameters:
- Compound: ${COMPOUND_NAME} (ID: ${COMPOUND_ID})
- Disease: ${DISEASE_CODE}
- Measures: Weight, Density
- Estimates: min, avg, max
- Lag: ${LAG_YEARS} years
- Models: M0, M1, M2, M3

Files Generated:
EOF

# Add file information
for file in ${RESULTS_DIR}/Results_*.csv; do
    if [ -f "$file" ]; then
        filename=$(basename "$file")
        linecount=$(($(wc -l < "$file") - 1))
        echo "- $filename: $linecount results" >> "$SUMMARY_FILE"
    fi
done

# Add comparison section
echo "" >> "$SUMMARY_FILE"
echo "Model Comparison Summary:" >> "$SUMMARY_FILE"
echo "------------------------" >> "$SUMMARY_FILE"

# Function to extract summary stats from CSV
extract_stats() {
    local file=$1
    local model_type=$2

    if [ -f "$file" ]; then
        echo "" >> "$SUMMARY_FILE"
        echo "${model_type} Models:" >> "$SUMMARY_FILE"

        # Count by dose-response type
        echo "  Dose-Response Types:" >> "$SUMMARY_FILE"
        awk -F',' 'NR>1 {count[$9]++} END {for (type in count) print "    " type ": " count[type]}' "$file" >> "$SUMMARY_FILE"

        # Average DIC by model
        echo "  Average DIC by Model:" >> "$SUMMARY_FILE"
        awk -F',' 'NR>1 && $17 != "NA" {sum[$8]+=$17; count[$8]++} END {for (model in sum) printf "    %s: %.2f\n", model, sum[model]/count[model]}' "$file" >> "$SUMMARY_FILE"

        # Count significant results
        echo "  Significant Results (p<0.05):" >> "$SUMMARY_FILE"
        awk -F',' 'NR>1 {if ($16 ~ /\*/) sig++; total++} END {printf "    %d/%d (%.1f%%)\n", sig, total, (sig/total)*100}' "$file" >> "$SUMMARY_FILE"

        # RR range for significant results
        echo "  RR Range (P90 vs P10, significant only):" >> "$SUMMARY_FILE"
        awk -F',' 'NR>1 && $16 ~ /\*/ && $13 != "NA" {
            if (min == "" || $13 < min) min = $13;
            if (max == "" || $13 > max) max = $13;
            sum += $13; count++
        } END {
            if (count > 0) {
                printf "    Min: %.4f, Max: %.4f, Mean: %.4f\n", min, max, sum/count
            } else {
                print "    No significant results"
            }
        }' "$file" >> "$SUMMARY_FILE"
    fi
}

# Extract stats for each analysis
extract_stats "${RESULTS_DIR}/Results_C81-C96_2,4-D_Linear_${TIMESTAMP}.csv" "Linear"
extract_stats "${RESULTS_DIR}/Results_C81-C96_2,4-D_Mixed_${TIMESTAMP}.csv" "Mixed"
extract_stats "${RESULTS_DIR}/Results_C81-C96_2,4-D_Nonlinear_${TIMESTAMP}.csv" "Non-linear"

# Create combined results file with best models
echo "" >> "$SUMMARY_FILE"
echo "Creating Combined Best Models File..." >> "$SUMMARY_FILE"

COMBINED_FILE="${RESULTS_DIR}/Combined_Best_Models.csv"

# Copy header from any results file
if [ -f "${RESULTS_DIR}/Results_C81-C96_2,4-D_Linear_${TIMESTAMP}.csv" ]; then
    head -n 1 "${RESULTS_DIR}/Results_C81-C96_2,4-D_Linear_${TIMESTAMP}.csv" > "$COMBINED_FILE"

    # For each combination, select model with lowest DIC
    echo "SELECT best models based on DIC..." >> "$SUMMARY_FILE"

    # This would require more complex processing - simplified version:
    # Combine all results and mark source
    for file in ${RESULTS_DIR}/Results_*.csv; do
        if [[ "$file" != *"Combined"* ]]; then
            tail -n +2 "$file" >> "${COMBINED_FILE}.tmp"
        fi
    done

    # Sort by measure, estimate, lag, and DIC to get best models
    sort -t',' -k5,5 -k6,6 -k7,7 -k17,17n "${COMBINED_FILE}.tmp" | \
    awk -F',' '!seen[$5","$6","$7]++' >> "$COMBINED_FILE"

    rm -f "${COMBINED_FILE}.tmp"

    echo "✓ Combined file created with best models" >> "$SUMMARY_FILE"
fi

# Display summary
echo ""
cat "$SUMMARY_FILE"

# Create visualization script
cat > "${RESULTS_DIR}/visualize_results.R" << 'EOF'
# R script to visualize 2,4-D analysis results
library(ggplot2)
library(dplyr)
library(readr)

# Read all result files
files <- list.files(pattern = "Results_.*\\.csv", full.names = TRUE)
results <- lapply(files, read_csv) %>% bind_rows()

# Create plots
# 1. DIC comparison
p1 <- ggplot(results %>% filter(!is.na(DIC)),
       aes(x = Model, y = DIC, color = Dose_Response_Type)) +
  geom_boxplot() +
  facet_grid(Measure ~ Estimate) +
  labs(title = "Model Fit Comparison (DIC)",
       subtitle = "2,4-D Analysis - Lower DIC indicates better fit") +
  theme_minimal()

# 2. RR comparison
p2 <- ggplot(results %>% filter(!is.na(RR_P90_vs_P10)),
       aes(x = Model, y = RR_P90_vs_P10, color = Dose_Response_Type)) +
  geom_point(position = position_dodge(0.3)) +
  geom_errorbar(aes(ymin = RR_P90_vs_P10_Lower, ymax = RR_P90_vs_P10_Upper),
                position = position_dodge(0.3), width = 0.2) +
  facet_grid(Measure ~ Estimate) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(title = "Relative Risk Comparison (P90 vs P10)",
       subtitle = "2,4-D Analysis - Error bars show 95% CI") +
  theme_minimal()

# Save plots
ggsave("DIC_comparison.png", p1, width = 10, height = 8)
ggsave("RR_comparison.png", p2, width = 10, height = 8)

print("Plots saved!")
EOF

echo ""
echo "📈 Visualization script created: ${RESULTS_DIR}/visualize_results.R"
echo "   Run in R to generate comparison plots"

# Final summary
echo ""
echo "============================================================"
echo "Analysis Complete!"
echo "============================================================"
echo "Results directory: ${RESULTS_DIR}"
echo ""
echo "Files created:"
ls -la "${RESULTS_DIR}/" | grep -E "\.(csv|txt)$"
echo ""
echo "Completed at: $(date)"
echo "============================================================"