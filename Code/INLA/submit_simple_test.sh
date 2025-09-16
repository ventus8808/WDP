#!/bin/bash

#================================================================================
# Slurm SBATCH Script for the Minimal INLA BYM2 Test
#================================================================================
# -- Job name
#SBATCH --job-name=INLA_BYM2_Test
#
# -- Output and error files
#SBATCH --output=slurm_logs/simple_test-%j.out
#SBATCH --error=slurm_logs/simple_test-%j.err
#
# -- Resource allocation
#SBATCH --partition=kshdtest
#SBATCH --gres=dcu:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:10:00 # 10 minutes is more than enough

#================================================================================
# Job Execution
#================================================================================
echo "========================================================"
echo "Starting Minimal INLA BYM2 Test"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Time: $(date)"
echo "========================================================"

# --- 1. Activate Conda Environment ---
echo "🚀 Activating Conda environment 'WDP'..."
source /public/home/acf4pijnzl/miniconda3/etc/profile.d/conda.sh
conda activate WDP
echo "   ✅ Conda environment activated."

# --- 2. Run the R Test Script ---
# We are already in the project root, so we specify the path to the script
echo "🚀 Running the R test script..."
Rscript Code/INLA/simple_bym2_test.R

echo "========================================================"
echo "Test script finished."
echo "Time: $(date)"
echo "========================================================"
