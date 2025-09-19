#!/bin/bash
# ========================
# WDP INLA Single Compound Analysis Script (Optimized Version)
# Robust script for running individual compound analysis
# ========================

#SBATCH --partition=kshctest
#SBATCH --job-name=WDP_INLA_Analysis
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=3:00:00
#SBATCH --output=%x-%j.log
#SBATCH --error=%x-%j.err
#SBATCH --exclude=a02r3n03,e08r3n02,e08r3n05,e08r3n09,e08r3n15,e10r3n00,e10r3n12,e10r4n11

# ========================
# Environment Setup
# ========================
echo "🚀 Starting WDP INLA Analysis"
echo "============================="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Job Name: ${SLURM_JOB_NAME}"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo ""

# Load required modules
echo "📦 Loading environment modules..."
module purge
module load compiler/dtk/23.10
module load rocm/5.3.3

# Activate Conda environment
echo "🐍 Activating Conda environment..."
source ~/miniconda3/etc/profile.d/conda.sh
conda activate INLA
echo "✅ Conda environment 'INLA' activated"
echo "R location: $(which Rscript)"
echo ""

# ========================
# Configuration
# ========================

# Project paths (updated for new structure)
PROJECT_ROOT="/public/home/acf4pijnzl/WDP"
SCRIPT_DIR="${PROJECT_ROOT}/Code/INLA"

# Default analysis parameters
DISEASE_CODE="${1:-C81-C96}"
COMPOUND_ID="${2:-1}"
MEASURE_TYPES="${3:-Weight,Density}"
ESTIMATE_TYPES="${4:-min}"
LAG_YEARS="${5:-5}"
MODEL_TYPES="${6:-M0,M1,M2,M3}"

echo "🔧 Analysis Configuration:"
echo "  Disease Code: ${DISEASE_CODE}"
echo "  Compound ID: ${COMPOUND_ID}"
echo "  Measure Types: ${MEASURE_TYPES}"
echo "  Estimate Types: ${ESTIMATE_TYPES}"
echo "  Lag Years: ${LAG_YEARS}"
echo "  Model Types: ${MODEL_TYPES}"
echo ""

# ========================
# Pre-flight Checks
# ========================

echo "🔍 Running pre-flight checks..."

# Check if project directory exists
if [ ! -d "${PROJECT_ROOT}" ]; then
    echo "❌ ERROR: Project directory not found: ${PROJECT_ROOT}"
    exit 1
fi

# Change to script directory
cd "${SCRIPT_DIR}" || {
    echo "❌ ERROR: Cannot change to script directory: ${SCRIPT_DIR}"
    exit 1
}

echo "📁 Working directory: $(pwd)"

# Check if main script exists
if [ ! -f "INLA_Main.R" ]; then
    echo "❌ ERROR: Main script not found: INLA_Main.R"
    exit 1
fi

# Check available disk space
echo "💾 Checking disk space..."
df -h . | tail -1 | awk '{print "Available space: " $4 " (Used: " $5 ")"}'

# Warn if less than 5GB available
AVAIL_KB=$(df . | tail -1 | awk '{print $4}')
if [ "$AVAIL_KB" -lt 5242880 ]; then  # 5GB in KB
    echo "⚠️  WARNING: Low disk space. INLA may fail if less than 5GB available."
fi

# Check dependencies
echo "🔍 Checking R package dependencies..."
if ! Rscript INLA_Dependencies/check_dependencies.R; then
    echo "❌ ERROR: Missing dependencies. Installing..."
    Rscript INLA_Dependencies/install_packages.R || {
        echo "❌ ERROR: Failed to install dependencies"
        exit 1
    }
fi

echo "✅ Pre-flight checks completed"
echo ""

# ========================
# Main Analysis Execution
# ========================

echo "🚀 Starting INLA analysis..."
echo "=============================="

# Construct the R command with all parameters
Rscript INLA_Main.R \
  --config "INLA_Config/analysis_config.yaml" \
  --disease-code "${DISEASE_CODE}" \
  --pesticide-category "compound:${COMPOUND_ID}" \
  --measure-type "${MEASURE_TYPES}" \
  --estimate-types "${ESTIMATE_TYPES}" \
  --lag-years "${LAG_YEARS}" \
  --model-types "${MODEL_TYPES}" \
  --verbose

# Capture exit status
STATUS=$?

echo ""
echo "=============================="
if [ $STATUS -eq 0 ]; then
    echo "✅ Analysis completed successfully!"
    echo "📊 Results saved to Result/INLA_Results/"
    echo "📋 Log files available in:"
    echo "  - Standard output: ${SLURM_JOB_NAME}-${SLURM_JOB_ID}.log"
    echo "  - Error output: ${SLURM_JOB_NAME}-${SLURM_JOB_ID}.err"
else
    echo "❌ Analysis failed with exit code: ${STATUS}"
    echo "🔍 Check error log: ${SLURM_JOB_NAME}-${SLURM_JOB_ID}.err"
    echo "💡 Common troubleshooting steps:"
    echo "  1. Check data file availability"
    echo "  2. Verify memory and time limits"
    echo "  3. Review parameter values"
fi

echo ""
echo "🕒 Job completed at: $(date)"
echo "⏱️  Total runtime: ${SECONDS} seconds"

exit ${STATUS}