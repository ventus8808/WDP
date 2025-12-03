#!/bin/bash
#SBATCH --partition=kshctest
#SBATCH --job-name=Ridge_C00_C97
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=Ridge_Test_C00_C97_%j.out
#SBATCH --error=Ridge_Test_C00_C97_%j.err

set -e

echo "========================================="
echo "Ridgeline Test - Simple Submit"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Start time: $(date)"
echo "========================================="

# Find project root
if [ -n "${SLURM_SUBMIT_DIR}" ]; then
  cd "$SLURM_SUBMIT_DIR"
fi

PROJECT_ROOT=$(pwd)
echo "Project root: $PROJECT_ROOT"

# Create output directories
mkdir -p Result/log
mkdir -p Result/Ridgeline
echo "Created output directories"

# Activate conda
ENV_NAME="brms"
echo "Activating conda environment: $ENV_NAME"

if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
  source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "/opt/anaconda3/etc/profile.d/conda.sh" ]; then
  source "/opt/anaconda3/etc/profile.d/conda.sh"
else
  echo "ERROR: Cannot find conda.sh"
  exit 1
fi

conda activate "$ENV_NAME" || {
  echo "ERROR: Failed to activate $ENV_NAME"
  exit 1
}

echo "Conda environment: $CONDA_DEFAULT_ENV"

# Load modules
module load devtoolset-8 2>/dev/null || echo "Warning: devtoolset-8 not loaded"

# Set environment
export TBB_CXX_TYPE=gcc
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo ""
echo "Environment variables:"
echo "  OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo "  MKL_NUM_THREADS=$MKL_NUM_THREADS"
echo ""

# Check script exists
SCRIPT="Code/brms/cmdstan_main_Ridgeline_Test.R"
if [ ! -f "$SCRIPT" ]; then
  echo "ERROR: Script not found: $SCRIPT"
  exit 1
fi
echo "Script: $SCRIPT"

# Check data exists
DATA="Data/Processed/df_EQI_AAMR_Triangulation/EQI_AAMR_Cluster_Climate.csv"
if [ ! -f "$DATA" ]; then
  echo "ERROR: Data not found: $DATA"
  exit 1
fi
echo "Data: $DATA"

echo ""
echo "========================================="
echo "Starting Ridgeline extraction..."
echo "========================================="
echo ""

# Run the script
Rscript "$SCRIPT" 2>&1 | tee "Result/log/Ridge_Test_${SLURM_JOB_ID}.log"

EXIT_CODE=${PIPESTATUS[0]}

echo ""
echo "========================================="
if [ $EXIT_CODE -eq 0 ]; then
  echo "SUCCESS!"
  echo "========================================="
  echo "Output: Result/Ridgeline/C00_C97_Ridge_Test.rds"

  if [ -f "Result/Ridgeline/C00_C97_Ridge_Test.rds" ]; then
    SIZE=$(du -h "Result/Ridgeline/C00_C97_Ridge_Test.rds" | cut -f1)
    echo "File size: $SIZE"
  fi

  echo ""
  echo "Next steps:"
  echo "  data <- readRDS('Result/Ridgeline/C00_C97_Ridge_Test.rds')"
  echo "  library(ggridges)"
  echo "  ggplot(data\$draws_long, aes(x=effect, y=quintile, fill=quintile)) +"
  echo "    geom_density_ridges(alpha=0.7)"
else
  echo "FAILED (exit code: $EXIT_CODE)"
  echo "Check log: Result/log/Ridge_Test_${SLURM_JOB_ID}.log"
fi
echo "========================================="
echo "End time: $(date)"

exit $EXIT_CODE
