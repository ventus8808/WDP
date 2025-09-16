#!/bin/bash

# Simple Docker run script for WDP INLA analysis
# Uses byminla-python-final:v4 image

# Default command (can be overridden by passing arguments)
DEFAULT_CMD="Rscript /project/Code/INLA/BYM_INLA_Production.R --help"
CMD="${@:-$DEFAULT_CMD}"

echo "Running in Docker: byminla-python-final:v4"
echo "Command: $CMD"
echo ""

# Run Docker container
docker run --rm \
  --platform linux/amd64 \
  -v "$(pwd)/../..:/project" \
  -w /project \
  byminla-python-final:v4 \
  bash -c "$CMD"