#!/bin/bash
# Local test wrapper for ridgeline extraction
# Run this to test the R script locally before submitting to SLURM

set -e

echo "========================================="
echo "Local Ridgeline Test"
echo "========================================="

# Find project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd -P)"

if [ ! -f "$PROJECT_ROOT/config.yaml" ]; then
  echo "ERROR: Cannot find config.yaml in $PROJECT_ROOT"
  exit 1
fi

cd "$PROJECT_ROOT"
echo "Project root: $PROJECT_ROOT"

# Create output directories
mkdir -p Result/log
mkdir -p Result/Ridgeline
echo "Created output directories"

# Check if R script exists
RUNNER="Code/brms/cmdstan_main_Ridgeline_Test.R"
if [ ! -f "$RUNNER" ]; then
  echo "ERROR: R script not found: $RUNNER"
  exit 1
fi

# Check if data exists
DATA_FILE="Data/Processed/df_EQI_AAMR_Triangulation/EQI_AAMR_Cluster_Climate.csv"
if [ ! -f "$DATA_FILE" ]; then
  echo "ERROR: Data file not found: $DATA_FILE"
  exit 1
fi

echo "Data file found: $DATA_FILE"
echo ""

# Activate conda environment
ENV_NAME="${ENV_NAME:-brms}"
echo "Activating conda environment: $ENV_NAME"

if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
  source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
  source "$HOME/anaconda3/etc/profile.d/conda.sh"
else
  echo "WARNING: Could not find conda initialization script"
fi

conda activate "$ENV_NAME" || {
  echo "ERROR: Failed to activate conda environment: $ENV_NAME"
  exit 1
}

echo "Conda environment: $CONDA_DEFAULT_ENV"
echo ""

# Run the R script
echo "========================================="
echo "Running Ridgeline Test R script..."
echo "========================================="
echo ""

Rscript "$RUNNER" 2>&1 | tee "Result/log/Ridge_Test_local_$(date +%Y%m%d_%H%M%S).log"

EXIT_CODE=${PIPESTATUS[0]}

echo ""
echo "========================================="
if [ $EXIT_CODE -eq 0 ]; then
  echo "✅ SUCCESS"
  echo "========================================="
  echo "Output: Result/Ridgeline/C00_C97_Ridge_Test.rds"

  if [ -f "Result/Ridgeline/C00_C97_Ridge_Test.rds" ]; then
    SIZE=$(du -h "Result/Ridgeline/C00_C97_Ridge_Test.rds" | cut -f1)
    echo "File size: $SIZE"
  fi
else
  echo "❌ FAILED (exit code: $EXIT_CODE)"
  echo "========================================="
  echo "Check log file in Result/log/"
fi

exit $EXIT_CODE
