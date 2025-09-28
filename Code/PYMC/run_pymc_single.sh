#!/bin/bash

# Local/cluster-agnostic runner for a single PyMC analysis
# Usage:
#   bash Code/PYMC/run_pymc_single.sh "C81-C96" "24D" "M3" "5" ["Weight"] [draws] [tune] [target_accept] [chains] [cores]
#
# Arguments:
#   1 disease (e.g., C81-C96)
#   2 compound (e.g., 24D)
#   3 model (M0|M1|M2|M3)
#   4 lag (e.g., 5)
#   5 measure (optional; default Weight)
#   6 draws (optional)
#   7 tune (optional)
#   8 target_accept (optional; e.g., 0.9)
#   9 chains (optional)
#   10 cores (optional)

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")"/../.. && pwd)"
PYMC_DIR="${PROJECT_ROOT}/Code/PYMC"
CONFIG_PATH="${PROJECT_ROOT}/config.yaml"

DISEASE=${1:?"Missing disease code"}
COMPOUND=${2:?"Missing compound"}
MODEL=${3:?"Missing model type"}
LAG=${4:?"Missing lag years"}
MEASURE=${5:-Weight}
DRAWS=${6:-}
TUNE=${7:-}
TARGET_ACCEPT=${8:-}
CHAINS=${9:-}
CORES=${10:-}

# Attempt to activate conda env 'pymc' if available
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
  source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "/opt/anaconda3/etc/profile.d/conda.sh" ]; then
  source "/opt/anaconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/opt/anaconda3/etc/profile.d/conda.sh" ]; then
  source "$HOME/opt/anaconda3/etc/profile.d/conda.sh"
fi

if command -v conda >/dev/null 2>&1 && conda info --envs | grep -q "^pymc\b"; then
  echo "Activating conda env: pymc"
  conda activate pymc || true
fi

cd "$PROJECT_ROOT"

CMD=(python "$PYMC_DIR/main.py" \
  --disease "$DISEASE" \
  --compound "$COMPOUND" \
  --model "$MODEL" \
  --lag "$LAG" \
  --measure "$MEASURE" \
  --sampling-mode production \
  --config-path "$CONFIG_PATH")

[ -n "$DRAWS" ] && CMD+=(--draws "$DRAWS")
[ -n "$TUNE" ] && CMD+=(--tune "$TUNE")
[ -n "$TARGET_ACCEPT" ] && CMD+=(--target-accept "$TARGET_ACCEPT")
[ -n "$CHAINS" ] && CMD+=(--chains "$CHAINS")
[ -n "$CORES" ] && CMD+=(--cores "$CORES")

echo "Running: ${CMD[*]}"
"${CMD[@]}"
