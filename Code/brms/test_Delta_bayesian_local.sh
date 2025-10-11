#!/bin/bash
# Local test script for Delta_bayesian.R
# Usage: bash Code/brms/test_Delta_bayesian_local.sh [CANCER_TYPE]
# Example: bash Code/brms/test_Delta_bayesian_local.sh C00_C97

set -eo pipefail

# Get cancer type from argument or use default
CANCER="${1:-C00_C97}"

echo "======================================================================"
echo "Testing Delta_bayesian.R locally"
echo "Cancer Type: $CANCER"
echo "======================================================================"

# Locate project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd -P)"
cd "$PROJECT_ROOT"

echo "Project root: $PROJECT_ROOT"

# Activate conda environment
ENV_NAME="${ENV_NAME:-brms}"
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "/opt/anaconda3/etc/profile.d/conda.sh" ]; then
    source "/opt/anaconda3/etc/profile.d/conda.sh"
else
    echo "ERROR: Cannot find conda initialization script"
    exit 1
fi

conda activate "$ENV_NAME" || { echo "ERROR: Failed to activate $ENV_NAME"; exit 1; }

echo "Conda environment: $CONDA_DEFAULT_ENV"
echo ""

# Run with test parameters (faster sampling for testing)
Rscript Code/brms/Delta_bayesian.R \
  --cancer-types "$CANCER" \
  --chains 2 \
  --iter 500 \
  --warmup 250 \
  --adapt-delta 0.95 \
  --max-treedepth 12 \
  --seed 1234

echo ""
echo "======================================================================"
echo "✅ Test complete!"
echo "Output file: Result/brms_delta/${CANCER}_delta.csv"
echo "======================================================================"
