#!/bin/bash

# Deploy and run WDP 2,4-D analysis on remote server
# Fixed version with correct paths and permissions

echo "============================================================"
echo "WDP Remote Server Deployment and Analysis (Fixed)"
echo "============================================================"
echo "Start time: $(date)"
echo ""

# Server configuration
SERVER_HOST="cancon.hpccube.com"
SERVER_PORT="65023"
SERVER_USER="acf4pijnzl"
SSH_KEY="acf4pijnzl_cancon.hpccube.com_RsaKeyExpireTime_2025-12-14_15-45-43.txt"

# Local paths
LOCAL_PROJECT_DIR="$(pwd)"
LOCAL_CODE_DIR="Code/INLA"
LOCAL_DATA_DIR="Data/Processed"

# Remote paths - FIXED to use correct home directory
REMOTE_HOME="/public/home/${SERVER_USER}"
REMOTE_PROJECT_DIR="${REMOTE_HOME}/WDP_Analysis"
REMOTE_RESULTS_DIR="${REMOTE_HOME}/WDP_Results"

# Analysis parameters
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
ANALYSIS_NAME="24D_Full_Analysis_${TIMESTAMP}"

# Check if SSH key exists
if [ ! -f "$SSH_KEY" ]; then
    echo "❌ Error: SSH key not found: $SSH_KEY"
    exit 1
fi

# Set proper permissions for SSH key
chmod 600 "$SSH_KEY"

# SSH command alias for convenience
SSH_CMD="ssh -o StrictHostKeyChecking=no -i $SSH_KEY -p $SERVER_PORT ${SERVER_USER}@${SERVER_HOST}"
SCP_CMD="scp -o StrictHostKeyChecking=no -i $SSH_KEY -P $SERVER_PORT"

echo "📡 Testing server connection..."
if $SSH_CMD "echo 'Connection successful'" 2>/dev/null; then
    echo "✓ Server connection established"
else
    echo "❌ Failed to connect to server"
    exit 1
fi

# Create remote directory structure
echo ""
echo "📁 Setting up remote directories..."
$SSH_CMD "mkdir -p ${REMOTE_PROJECT_DIR}/{Code/INLA/{utils,config},Data/Processed/{CDC,Pesticide,PCA,Socioeconomic},Results,Config} && \
          mkdir -p ${REMOTE_RESULTS_DIR}/${ANALYSIS_NAME} && \
          echo '✓ Remote directories created successfully'"

# Transfer code files
echo ""
echo "📤 Transferring analysis code..."

# Create tar archive of code files
tar -czf /tmp/wdp_code.tar.gz -C ${LOCAL_CODE_DIR} . 2>/dev/null

# Transfer and extract code
$SCP_CMD /tmp/wdp_code.tar.gz ${SERVER_USER}@${SERVER_HOST}:/tmp/
$SSH_CMD "cd ${REMOTE_PROJECT_DIR}/Code/INLA && tar -xzf /tmp/wdp_code.tar.gz && rm /tmp/wdp_code.tar.gz && echo '✓ Code files transferred and extracted'"

# Clean up local archive
rm /tmp/wdp_code.tar.gz

# Check if data needs to be transferred
echo ""
echo "📊 Checking data files on server..."
DATA_EXISTS=$($SSH_CMD "if [ -d '${REMOTE_PROJECT_DIR}/Data/Processed/CDC' ] && [ \"\$(ls -A ${REMOTE_PROJECT_DIR}/Data/Processed/CDC 2>/dev/null)\" ]; then echo 'YES'; else echo 'NO'; fi")

if [ "$DATA_EXISTS" = "NO" ]; then
    echo "📤 Transferring data files (this may take a while)..."

    # Transfer data files by directory to avoid large archives
    for subdir in CDC Pesticide PCA Socioeconomic; do
        if [ -d "${LOCAL_DATA_DIR}/${subdir}" ]; then
            echo "   Transferring ${subdir} data..."
            tar -czf /tmp/wdp_${subdir}.tar.gz -C ${LOCAL_DATA_DIR}/${subdir} .
            $SCP_CMD /tmp/wdp_${subdir}.tar.gz ${SERVER_USER}@${SERVER_HOST}:/tmp/
            $SSH_CMD "cd ${REMOTE_PROJECT_DIR}/Data/Processed/${subdir} && tar -xzf /tmp/wdp_${subdir}.tar.gz && rm /tmp/wdp_${subdir}.tar.gz"
            rm /tmp/wdp_${subdir}.tar.gz
        fi
    done
    echo "✓ Data files transferred"
else
    echo "✓ Data files already exist on server"
fi

# Create analysis script on server
echo ""
echo "📝 Creating remote analysis script..."

# Use a here document to create the script
$SSH_CMD "cat > ${REMOTE_PROJECT_DIR}/run_24d_analysis_server.sh" << 'EOF'
#!/bin/bash

# Server-side analysis script for 2,4-D
# This runs directly on the HPC server

echo "============================================================"
echo "WDP 2,4-D Analysis on HPC Server"
echo "============================================================"
echo "Server: $(hostname)"
echo "Start time: $(date)"
echo ""

# Set analysis parameters
COMPOUND_ID="2"
COMPOUND_NAME="2,4-D"
DISEASE_CODE="C81-C96"
MEASURE_TYPES="Weight,Density"
ESTIMATE_TYPES="min,avg,max"
LAG_YEARS="5"
MODEL_TYPES="M0,M1,M2,M3"

# Results directory (passed as argument)
RESULTS_DIR="$1"
PROJECT_DIR="$2"

# Change to project directory
cd "$PROJECT_DIR"

# Load R module if available
if command -v module &> /dev/null; then
    echo "Loading R module..."
    module load R/4.3.0 2>/dev/null || module load R 2>/dev/null || echo "Using system R"
fi

# Check R installation
echo "R version:"
R --version | head -n 1

# Install required packages if not available
echo ""
echo "📦 Checking R packages..."
R --quiet --slave << 'RSCRIPT'
# Function to check and install packages
check_packages <- function() {
    # Set CRAN mirror
    options(repos = c(CRAN = "https://cran.rstudio.com/"))

    required_packages <- c("dplyr", "readr", "yaml", "argparse", "progress")

    # First check and install INLA separately
    if (!requireNamespace("INLA", quietly = TRUE)) {
        cat("Installing INLA package...\n")
        install.packages("INLA", repos = c(getOption("repos"),
                         INLA = "https://inla.r-inla-download.org/R/stable"))
    }

    # Check other packages
    missing <- required_packages[!required_packages %in% installed.packages()[,"Package"]]

    if (length(missing) > 0) {
        cat("Installing missing packages:", paste(missing, collapse=", "), "\n")
        install.packages(missing)
    }

    # Verify all packages can be loaded
    all_packages <- c("INLA", required_packages)
    for (pkg in all_packages) {
        if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
            stop(paste("Failed to load package:", pkg))
        }
    }

    cat("✓ All required packages are available\n")
}

# Run package check
tryCatch({
    check_packages()
}, error = function(e) {
    cat("Error installing packages:", e$message, "\n")
    cat("Proceeding anyway...\n")
})
RSCRIPT

# Create config files for different model types
CONFIG_DIR="${RESULTS_DIR}/configs"
mkdir -p "$CONFIG_DIR"

# Copy and modify configs
cp Code/INLA/config/analysis_config.yaml "$CONFIG_DIR/config_linear.yaml"

# Fix paths in config to use /project instead of absolute paths
sed -i 's|base_dir: "/project/Data/Processed"|base_dir: "Data/Processed"|g' "$CONFIG_DIR/config_linear.yaml"
sed -i 's|base_dir: "/project/Result/Filter"|base_dir: "Results"|g' "$CONFIG_DIR/config_linear.yaml"

# Create non-linear config
sed 's/enabled: false # Set to true to enable non-linear models/enabled: true # Set to true to enable non-linear models/' \
    "$CONFIG_DIR/config_linear.yaml" > "$CONFIG_DIR/config_nonlinear_all.yaml"
sed -i 's/model_types: \["M2", "M3"\]/model_types: ["M0", "M1", "M2", "M3"]/' "$CONFIG_DIR/config_nonlinear_all.yaml"

# Create mixed config
sed 's/enabled: false # Set to true to enable non-linear models/enabled: true # Set to true to enable non-linear models/' \
    "$CONFIG_DIR/config_linear.yaml" > "$CONFIG_DIR/config_nonlinear_m2m3.yaml"

# Function to run analysis
run_analysis() {
    local config_file=$1
    local output_name=$2
    local description=$3

    echo ""
    echo "🚀 Running: ${description}"
    echo "   Config: ${config_file}"
    echo "   Output: ${output_name}"

    # Run R script
    Rscript Code/INLA/BYM_INLA_Production.R \
        --config "${CONFIG_DIR}/${config_file}" \
        --pesticide-category "compound:${COMPOUND_ID}" \
        --measure-type "${MEASURE_TYPES}" \
        --estimate-types "${ESTIMATE_TYPES}" \
        --lag-years "${LAG_YEARS}" \
        --model-types "${MODEL_TYPES}" \
        --disease-code "${DISEASE_CODE}" \
        --output-file "${RESULTS_DIR}/${output_name}" \
        --verbose 2>&1 | tee "${RESULTS_DIR}/${output_name}.log"

    # Check results
    if [ -f "${RESULTS_DIR}/${output_name}" ]; then
        local num_results=$(wc -l < "${RESULTS_DIR}/${output_name}")
        echo "   ✓ Completed: $((num_results - 1)) results"
    else
        echo "   ❌ Failed - check log file"
    fi
}

# Run analyses
echo ""
echo "🔄 Starting Analysis Runs"
echo "========================"

# Run all three model types
run_analysis "config_linear.yaml" \
    "Results_24D_Linear.csv" \
    "Linear dose-response models"

run_analysis "config_nonlinear_m2m3.yaml" \
    "Results_24D_Mixed.csv" \
    "Mixed models (Linear M0,M1 + Non-linear M2,M3)"

run_analysis "config_nonlinear_all.yaml" \
    "Results_24D_Nonlinear.csv" \
    "Non-linear dose-response models"

# Generate summary
echo ""
echo "📊 Generating Summary"
echo "==================="

# Create summary file
cat > "${RESULTS_DIR}/Analysis_Summary.txt" << SUMMARY
WDP 2,4-D Analysis Summary
==========================
Server: $(hostname)
Completed: $(date)

Parameters:
- Compound: ${COMPOUND_NAME} (ID: ${COMPOUND_ID})
- Disease: ${DISEASE_CODE}
- Measures: ${MEASURE_TYPES}
- Estimates: ${ESTIMATE_TYPES}
- Lag: ${LAG_YEARS} years
- Models: ${MODEL_TYPES}

Results Files:
SUMMARY

# Add file stats
for file in ${RESULTS_DIR}/Results_*.csv; do
    if [ -f "$file" ]; then
        filename=$(basename "$file")
        lines=$(($(wc -l < "$file") - 1))
        size=$(du -h "$file" | cut -f1)
        echo "- $filename: $lines results ($size)" >> "${RESULTS_DIR}/Analysis_Summary.txt"
    fi
done

# Add basic statistics
echo "" >> "${RESULTS_DIR}/Analysis_Summary.txt"
echo "Model Performance Summary:" >> "${RESULTS_DIR}/Analysis_Summary.txt"
echo "=========================" >> "${RESULTS_DIR}/Analysis_Summary.txt"

for file in ${RESULTS_DIR}/Results_*.csv; do
    if [ -f "$file" ] && [ $(wc -l < "$file") -gt 1 ]; then
        filename=$(basename "$file" .csv)
        echo "" >> "${RESULTS_DIR}/Analysis_Summary.txt"
        echo "$filename:" >> "${RESULTS_DIR}/Analysis_Summary.txt"

        # Count dose-response types
        awk -F',' 'NR>1 {count[$9]++} END {for (type in count) print "  " type ": " count[type] " models"}' "$file" >> "${RESULTS_DIR}/Analysis_Summary.txt"

        # Count significant results
        awk -F',' 'NR>1 {total++; if ($16 ~ /\*/) sig++} END {printf "  Significant results: %d/%d (%.1f%%)\n", sig, total, (sig/total)*100}' "$file" >> "${RESULTS_DIR}/Analysis_Summary.txt"
    fi
done

echo ""
echo "✓ Analysis complete!"
echo "Results saved to: ${RESULTS_DIR}"
echo "Completed at: $(date)"
echo "============================================================"
EOF

# Make script executable
$SSH_CMD "chmod +x ${REMOTE_PROJECT_DIR}/run_24d_analysis_server.sh"
echo "✓ Remote analysis script created"

# Execute analysis on server
echo ""
echo "🚀 Starting analysis on server..."
echo "   This may take 30-60 minutes depending on server load"
echo ""

# Run with nohup to prevent disconnection issues
$SSH_CMD "cd ${REMOTE_PROJECT_DIR} && \
          nohup bash run_24d_analysis_server.sh '${REMOTE_RESULTS_DIR}/${ANALYSIS_NAME}' '${REMOTE_PROJECT_DIR}' \
          > '${REMOTE_RESULTS_DIR}/${ANALYSIS_NAME}/server_execution.log' 2>&1 & \
          echo \$! > '${REMOTE_RESULTS_DIR}/${ANALYSIS_NAME}/analysis.pid' && \
          echo '✓ Analysis started with PID: '\$(cat '${REMOTE_RESULTS_DIR}/${ANALYSIS_NAME}/analysis.pid')"

# Wait a moment and show initial output
sleep 5
echo ""
echo "Initial output:"
echo "---------------"
$SSH_CMD "head -n 30 '${REMOTE_RESULTS_DIR}/${ANALYSIS_NAME}/server_execution.log' 2>/dev/null || echo 'Waiting for log file...'"

# Create monitoring script
cat > monitor_analysis.sh << MONITOR
#!/bin/bash
# Monitor script to check analysis progress

SERVER_HOST="cancon.hpccube.com"
SERVER_PORT="65023"
SERVER_USER="acf4pijnzl"
SSH_KEY="acf4pijnzl_cancon.hpccube.com_RsaKeyExpireTime_2025-12-14_15-45-43.txt"
ANALYSIS_NAME="${ANALYSIS_NAME}"

SSH_CMD="ssh -o StrictHostKeyChecking=no -i \$SSH_KEY -p \$SERVER_PORT \${SERVER_USER}@\${SERVER_HOST}"
REMOTE_HOME="/public/home/\${SERVER_USER}"
REMOTE_RESULTS="\${REMOTE_HOME}/WDP_Results/\${ANALYSIS_NAME}"

echo "Monitoring analysis: \$ANALYSIS_NAME"
echo "Time: \$(date)"
echo ""

# Check if process is still running
\$SSH_CMD << EOF
    if [ -f "\${REMOTE_RESULTS}/analysis.pid" ]; then
        PID=\\\$(cat "\${REMOTE_RESULTS}/analysis.pid")
        if ps -p \\\$PID > /dev/null 2>&1; then
            echo "✓ Analysis is still running (PID: \\\$PID)"
            echo "  Running for: \\\$(ps -o etime= -p \\\$PID | xargs)"
        else
            echo "✓ Analysis has completed"
        fi
    else
        echo "⚠️  PID file not found"
    fi

    echo ""
    echo "Recent output (last 30 lines):"
    echo "==============================="
    tail -n 30 "\${REMOTE_RESULTS}/server_execution.log" 2>/dev/null || echo "Log file not found yet"

    echo ""
    echo "Result files:"
    echo "============="
    ls -lah "\${REMOTE_RESULTS}"/Results_*.csv 2>/dev/null || echo "No result files yet"

    echo ""
    echo "Log files:"
    echo "=========="
    ls -lah "\${REMOTE_RESULTS}"/*.log 2>/dev/null | grep -v server_execution || echo "No log files yet"
EOF
MONITOR

chmod +x monitor_analysis.sh

# Create download script
cat > download_results.sh << DOWNLOAD
#!/bin/bash
# Download results from server

SERVER_HOST="cancon.hpccube.com"
SERVER_PORT="65023"
SERVER_USER="acf4pijnzl"
SSH_KEY="acf4pijnzl_cancon.hpccube.com_RsaKeyExpireTime_2025-12-14_15-45-43.txt"
ANALYSIS_NAME="${ANALYSIS_NAME}"

REMOTE_HOME="/public/home/\${SERVER_USER}"
REMOTE_RESULTS="\${REMOTE_HOME}/WDP_Results/\${ANALYSIS_NAME}"
LOCAL_RESULTS="Results/\${ANALYSIS_NAME}"

echo "Downloading results for: \$ANALYSIS_NAME"
mkdir -p "\$LOCAL_RESULTS"

# Use rsync to download all results
rsync -avz --progress \
    -e "ssh -o StrictHostKeyChecking=no -i \$SSH_KEY -p \$SERVER_PORT" \
    \${SERVER_USER}@\${SERVER_HOST}:\${REMOTE_RESULTS}/ \${LOCAL_RESULTS}/

echo ""
echo "✓ Results downloaded to: \$LOCAL_RESULTS"
echo ""
echo "Files downloaded:"
ls -lah "\$LOCAL_RESULTS"
echo ""
echo "Summary:"
if [ -f "\$LOCAL_RESULTS/Analysis_Summary.txt" ]; then
    cat "\$LOCAL_RESULTS/Analysis_Summary.txt"
fi
DOWNLOAD

chmod +x download_results.sh

# Final instructions
echo ""
echo "============================================================"
echo "✓ Deployment Complete!"
echo "============================================================"
echo ""
echo "Analysis Name: ${ANALYSIS_NAME}"
echo ""
echo "The analysis is now running on the server."
echo "Expected completion time: 30-60 minutes"
echo ""
echo "📋 Quick Commands:"
echo "  Monitor progress:  ./monitor_analysis.sh"
echo "  Download results:  ./download_results.sh"
echo ""
echo "📍 Server Locations:"
echo "  Project: ${REMOTE_PROJECT_DIR}"
echo "  Results: ${REMOTE_RESULTS_DIR}/${ANALYSIS_NAME}"
echo ""
echo "The analysis will run all combinations:"
echo "  - 2 Measures (Weight, Density)"
echo "  - 3 Estimates (min, avg, max)"
echo "  - 4 Models (M0, M1, M2, M3)"
echo "  - 3 Dose-response types (Linear, Mixed, Non-linear)"
echo "  = Total: 72 model runs"
echo ""
echo "============================================================"
