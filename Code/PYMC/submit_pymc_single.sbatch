#!/bin/bash
#SBATCH --job-name=pymc_single
#SBATCH --output=logs/pymc_single_%j.out
#SBATCH --error=logs/pymc_single_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

# SBATCH usage:
# sbatch Code/PYMC/submit_pymc_single.sbatch C81-C96 24D M3 5 Weight 4000 2000 0.9 4 8

set -euo pipefail

DISEASE=${1:?Missing disease}
COMPOUND=${2:?Missing compound}
MODEL=${3:?Missing model}
LAG=${4:?Missing lag}
MEASURE=${5:-Weight}
DRAWS=${6:-4000}
TUNE=${7:-2000}
TARGET_ACCEPT=${8:-0.9}
CHAINS=${9:-4}
CORES=${10:-$SLURM_CPUS_PER_TASK}

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")"/../.. && pwd)"
RUNNER="${PROJECT_ROOT}/Code/PYMC/run_pymc_single.sh"

mkdir -p "${PROJECT_ROOT}/logs"

# Ensure conda is available on cluster nodes
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
  source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "/opt/anaconda3/etc/profile.d/conda.sh" ]; then
  source "/opt/anaconda3/etc/profile.d/conda.sh"
fi

# Activate environment if present
if command -v conda >/dev/null 2>&1 && conda info --envs | grep -q "^pymc\b"; then
  conda activate pymc || true
fi

bash "$RUNNER" "$DISEASE" "$COMPOUND" "$MODEL" "$LAG" "$MEASURE" "$DRAWS" "$TUNE" "$TARGET_ACCEPT" "$CHAINS" "$CORES"
