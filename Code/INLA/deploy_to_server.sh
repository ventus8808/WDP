#!/bin/bash

# Deploy and run WDP 2,4-D analysis on remote server
# This script transfers files and executes analysis on HPC server

echo "============================================================"
echo "WDP Remote Server Deployment and Analysis"
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

# Remote paths
REMOTE_HOME="/home/${SERVER_USER}"
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
SSH_CMD="ssh -i $SSH_KEY -p $SERVER_PORT ${SERVER_USER}@${SERVER_HOST}"
SCP_CMD="scp -i $SSH_KEY -P $SERVER_PORT"
RSYNC_CMD="rsync -avz -e 'ssh -i $SSH_KEY -p $SERVER_PORT'"

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
$SSH_CMD << EOF
    mkdir -p ${REMOTE_PROJECT_DIR}/{Code/INLA/utils,Data/Processed,Results,Config}
    mkdir -p ${REMOTE_RESULTS_DIR}/${ANALYSIS_NAME}
    echo "✓ Remote directories created"
EOF

# Transfer code files
echo ""
echo "📤 Transferring analysis code..."

# Transfer main scripts
$RSYNC_CMD \
    ${LOCAL_CODE_DIR}/*.R \
    ${LOCAL_CODE_DIR}/*.sh \
    ${SERVER_USER}@${SERVER_HOST}:${REMOTE_PROJECT_DIR}/Code/INLA/

# Transfer utility modules
$RSYNC_CMD \
    ${LOCAL_CODE_DIR}/utils/*.R \
    ${SERVER_USER}@${SERVER_HOST}:${REMOTE_PROJECT_DIR}/Code/INLA/utils/

# Transfer configuration
$RSYNC_CMD \
    ${LOCAL_CODE_DIR}/config/analysis_config.yaml \
    ${SERVER_USER}@${SERVER_HOST}:${REMOTE_PROJECT_DIR}/Code/INLA/config/

echo "✓ Code files transferred"

# Check if data needs to be transferred
echo ""
echo "📊 Checking data files on server..."
$SSH_CMD << EOF
    if [ -d "${REMOTE_PROJECT_DIR}/Data/Processed" ] && [ "\$(ls -A ${REMOTE_PROJECT_DIR}/Data/Processed 2>/dev/null)" ]; then
        echo "✓ Data files already exist on server"
        exit 0
    else
        echo "⚠️  Data files need to be transferred"
        exit 1
    fi
EOF

if [ $? -ne 0 ]; then
    echo "📤 Transferring data files (this may take a while)..."

    # Create a compressed archive of data files
    echo "   Creating data archive..."
    tar -czf /tmp/wdp_data.tar.gz -C ${LOCAL_DATA_DIR} .

    # Transfer the archive
    echo "   Transferring archive..."
    $SCP_CMD /tmp/wdp_data.tar.gz ${SERVER_USER}@${SERVER_HOST}:/tmp/

    # Extract on server
    echo "   Extracting data on server..."
    $SSH_CMD << EOF
        cd ${REMOTE_PROJECT_DIR}/Data/Processed
        tar -xzf /tmp/wdp_data.tar.gz
        rm /tmp/wdp_data.tar.gz
        echo "✓ Data files extracted"
EOF

    # Clean up local archive
    rm /tmp/wdp_data.tar.gz
fi

# Create analysis script on server
echo ""
echo "📝 Creating remote analysis script..."

$SSH_CMD << 'REMOTE_SCRIPT'
cat > ${REMOTE_PROJECT_DIR}/run_24d_analysis_server.sh << 'EOF'
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

# Load R module if available
if command -v module &> /dev/null; then
    echo "Loading R module..."
    module load R/4.3.0 || module load R || echo "R module not found, using system R"
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

# Create config files for different model types
CONFIG_DIR="${RESULTS_DIR}/configs"
mkdir -p "$CONFIG_DIR"

# Copy and modify configs
cp Code/INLA/config/analysis_config.yaml "$CONFIG_DIR/config_linear.yaml"

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

# Linear models
run_analysis "config_linear.yaml" \
    "Results_24D_Linear.csv" \
    "Linear dose-response models"

# Mixed models
run_analysis "config_nonlinear_m2m3.yaml" \
    "Results_24D_Mixed.csv" \
    "Mixed models (Linear M0,M1 + Non-linear M2,M3)"

# Full non-linear
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
        echo "- $filename: $lines results" >> "${RESULTS_DIR}/Analysis_Summary.txt"
    fi
done

echo ""
echo "✓ Analysis complete!"
echo "Results saved to: ${RESULTS_DIR}"
echo "Completed at: $(date)"
echo "============================================================"
EOF

# Make script executable
chmod +x ${REMOTE_PROJECT_DIR}/run_24d_analysis_server.sh
echo "✓ Remote analysis script created"
REMOTE_SCRIPT

# Execute analysis on server
echo ""
echo "🚀 Starting analysis on server..."
echo "   This may take 30-60 minutes depending on server load"
echo ""

# Run with nohup to prevent disconnection issues
$SSH_CMD << EOF
    cd ${REMOTE_PROJECT_DIR}
    nohup bash run_24d_analysis_server.sh "${REMOTE_RESULTS_DIR}/${ANALYSIS_NAME}" \
        > "${REMOTE_RESULTS_DIR}/${ANALYSIS_NAME}/server_execution.log" 2>&1 &

    # Get the process ID
    PID=\$!
    echo "✓ Analysis started with PID: \$PID"
    echo \$PID > "${REMOTE_RESULTS_DIR}/${ANALYSIS_NAME}/analysis.pid"

    # Show initial log output
    sleep 5
    echo ""
    echo "Initial output:"
    echo "---------------"
    head -n 20 "${REMOTE_RESULTS_DIR}/${ANALYSIS_NAME}/server_execution.log"
EOF

# Create monitoring script
cat > monitor_analysis.sh << 'MONITOR'
#!/bin/bash
# Monitor script to check analysis progress

SERVER_HOST="cancon.hpccube.com"
SERVER_PORT="65023"
SERVER_USER="acf4pijnzl"
SSH_KEY="acf4pijnzl_cancon.hpccube.com_RsaKeyExpireTime_2025-12-14_15-45-43.txt"
ANALYSIS_NAME="$1"

if [ -z "$ANALYSIS_NAME" ]; then
    echo "Usage: ./monitor_analysis.sh ANALYSIS_NAME"
    exit 1
fi

SSH_CMD="ssh -i $SSH_KEY -p $SERVER_PORT ${SERVER_USER}@${SERVER_HOST}"
REMOTE_RESULTS="/home/${SERVER_USER}/WDP_Results/${ANALYSIS_NAME}"

echo "Monitoring analysis: $ANALYSIS_NAME"
echo ""

# Check if process is still running
$SSH_CMD << EOF
    if [ -f "${REMOTE_RESULTS}/analysis.pid" ]; then
        PID=\$(cat "${REMOTE_RESULTS}/analysis.pid")
        if ps -p \$PID > /dev/null; then
            echo "✓ Analysis is still running (PID: \$PID)"
        else
            echo "✓ Analysis has completed"
        fi
    fi

    echo ""
    echo "Recent output:"
    echo "-------------"
    tail -n 20 "${REMOTE_RESULTS}/server_execution.log" 2>/dev/null || echo "Log file not found yet"

    echo ""
    echo "Result files:"
    echo "------------"
    ls -la "${REMOTE_RESULTS}"/Results_*.csv 2>/dev/null || echo "No result files yet"
EOF
MONITOR

chmod +x monitor_analysis.sh

# Create download script
cat > download_results.sh << 'DOWNLOAD'
#!/bin/bash
# Download results from server

SERVER_HOST="cancon.hpccube.com"
SERVER_PORT="65023"
SERVER_USER="acf4pijnzl"
SSH_KEY="acf4pijnzl_cancon.hpccube.com_RsaKeyExpireTime_2025-12-14_15-45-43.txt"
ANALYSIS_NAME="$1"

if [ -z "$ANALYSIS_NAME" ]; then
    echo "Usage: ./download_results.sh ANALYSIS_NAME"
    exit 1
fi

RSYNC_CMD="rsync -avz -e 'ssh -i $SSH_KEY -p $SERVER_PORT'"
REMOTE_RESULTS="/home/${SERVER_USER}/WDP_Results/${ANALYSIS_NAME}"
LOCAL_RESULTS="Results/${ANALYSIS_NAME}"

echo "Downloading results for: $ANALYSIS_NAME"
mkdir -p "$LOCAL_RESULTS"

$RSYNC_CMD ${SERVER_USER}@${SERVER_HOST}:${REMOTE_RESULTS}/ ${LOCAL_RESULTS}/

echo ""
echo "✓ Results downloaded to: $LOCAL_RESULTS"
ls -la "$LOCAL_RESULTS"
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
echo "This will take approximately 30-60 minutes to complete."
echo ""
echo "To monitor progress:"
echo "  ./monitor_analysis.sh ${ANALYSIS_NAME}"
echo ""
echo "To download results when complete:"
echo "  ./download_results.sh ${ANALYSIS_NAME}"
echo ""
echo "Server paths:"
echo "  Project: ${REMOTE_PROJECT_DIR}"
echo "  Results: ${REMOTE_RESULTS_DIR}/${ANALYSIS_NAME}"
echo ""
echo "============================================================"
