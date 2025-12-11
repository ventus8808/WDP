#!/bin/bash
#SBATCH --partition=kshctest
#SBATCH --job-name=main_ridgeline
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=1-00:00:00
#SBATCH --output=main_ridgeline_%A_%a.out
#SBATCH --error=main_ridgeline_%A_%a.err

# SLURM array job for ridgeline posterior extraction
# Array tasks: 3 jobs (one per lag year, each runs 2 models)
#
# Task mapping:
#   1 → Lag5  (runs Overall + MultiDomain)
#   2 → Lag10 (runs Overall + MultiDomain)
#   3 → Lag15 (runs Overall + MultiDomain)
#
# Usage: sbatch submit_ridgeline.sh

set -e

PROJECT_ROOT="${PROJECT_ROOT:-/home/$USER/WDP}"
cd "$PROJECT_ROOT"

mkdir -p logs

RSCRIPT="Code/brms/cmdstan_main_ridgeline.R"

if [ ! -f "$RSCRIPT" ]; then
    echo "ERROR: R script not found: $RSCRIPT"
    exit 1
fi

# Load R environment
if [ -f "$HOME/miniforge3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniforge3/etc/profile.d/conda.sh"
    conda activate wdp
elif [ -f "$HOME/mambaforge/etc/profile.d/mamba.sh" ]; then
    source "$HOME/mambaforge/etc/profile.d/mamba.sh"
    mamba activate wdp
fi

TASK_ID=${SLURM_ARRAY_TASK_ID}
SCENARIO=$TASK_ID

case $TASK_ID in
    1) DESC="Lag5" ;;
    2) DESC="Lag10" ;;
    3) DESC="Lag15" ;;
    *)
        echo "ERROR: Invalid SLURM_ARRAY_TASK_ID: $TASK_ID"
        exit 1
        ;;
esac

echo "========================================"
echo "Ridgeline Production - $DESC"
echo "========================================"
echo "Job ID:    $SLURM_JOB_ID"
echo "Task:      $TASK_ID / 3"
echo "Scenario:  $SCENARIO"
echo "Node:      $SLURMD_NODENAME"
echo "Started:   $(date)"
echo "========================================"

CANCER="C00_C97"
CHAINS=4
ITER=2000
WARMUP=1000
ADAPT_DELTA=0.95
MAX_TREEDEPTH=12
SEED=1234

# Run Overall model
echo ""
echo "Running Overall EQI model..."
Rscript "$RSCRIPT" \
    --scenario "$SCENARIO" \
    --model "overall" \
    --cancer "$CANCER" \
    --chains "$CHAINS" \
    --iter "$ITER" \
    --warmup "$WARMUP" \
    --adapt-delta "$ADAPT_DELTA" \
    --max-treedepth "$MAX_TREEDEPTH" \
    --seed "$SEED"

EXIT1=$?

if [ $EXIT1 -ne 0 ]; then
    echo "ERROR: Overall model failed"
    exit $EXIT1
fi

# Run Multi-domain model
echo ""
echo "Running Multi-domain model..."
Rscript "$RSCRIPT" \
    --scenario "$SCENARIO" \
    --model "multi" \
    --cancer "$CANCER" \
    --chains "$CHAINS" \
    --iter "$ITER" \
    --warmup "$WARMUP" \
    --adapt-delta "$ADAPT_DELTA" \
    --max-treedepth "$MAX_TREEDEPTH" \
    --seed "$SEED"

EXIT2=$?

if [ $EXIT2 -ne 0 ]; then
    echo "ERROR: Multi-domain model failed"
    exit $EXIT2
fi

echo ""
echo "========================================"
echo "✓ Task $TASK_ID completed successfully"
echo "  Lag: $DESC"
echo "  Models: Overall + MultiDomain"
echo "Finished: $(date)"
echo "========================================"

exit 0
