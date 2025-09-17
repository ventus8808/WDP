#!/bin/bash
# Smoke test for the run_compound wrapper (dry-run)
set -euo pipefail

ROOT=$(dirname "$0")/..
ROOT=$(cd "$ROOT" && pwd)

echo "Running smoke checks..."

python3 tools/run_compound.py --jobs examples/compound_jobs.yaml --dry-run

echo "Checking that run and submit scripts exist..."
for f in Code/INLA/run_single_compound_test.sh Code/INLA/submit_single_compound_test.sh; do
  if [ -f "$f" ]; then
    echo "  ✓ $f"
  else
    echo "  ✗ $f (MISSING)"
    exit 1
  fi
done

echo "Smoke run complete. To run locally:"
echo "  python3 tools/run_compound.py --compound-id 5 --compound-name Abamectin --disease C81-C96 --local"
echo "To submit (server):"
echo "  python3 tools/run_compound.py --jobs examples/compound_jobs.yaml --submit"
