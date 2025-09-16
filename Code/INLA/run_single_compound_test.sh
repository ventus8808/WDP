#!/bin/bash

# WDP Single Compound Test Script (Server Version - No Docker)
# Tests all model combinations for a single compound
# Author: WDP Analysis Team
# Date: 2024

echo "============================================================"
echo "WDP Single Compound Analysis Test (Server Version)"
echo "============================================================"
echo "Start time: $(date)"
echo ""

# Set default parameters (can be overridden via command line)
COMPOUND_ID="${1:-2}"           # Default: 2,4-D (compound ID 2)
COMPOUND_NAME="${2:-2,4-D}"     # Default: 2,4-D name
DISEASE_CODE="${3:-C81-C96}"    # Default: Lymphoid and Hematopoietic
MEASURE_TYPES="${4:-Weight,Density}"     # Both weight and density
ESTIMATE_TYPES="${5:-min,avg,max}"       # All estimates
LAG_YEARS="${6:-5}"             # 5-year lag
MODEL_TYPES="${7:-M0,M1,M2,M3}" # All models

# Create timestamp for output
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

echo "📋 Analysis Parameters:"
echo "----------------------"
echo "  Compound: ${COMPOUND_NAME} (ID: ${COMPOUND_ID})"
echo "  Disease: ${DISEASE_CODE}"
echo "  Measures: ${MEASURE_TYPES}"
echo "  Estimates: ${ESTIMATE_TYPES}"
echo "  Lag: ${LAG_YEARS} years"
echo "  Models: ${MODEL_TYPES}"
echo "  Timestamp: ${TIMESTAMP}"
echo ""

# Set working directory (adjust this path according to your server setup)
WORK_DIR="$(pwd)"
echo "Working directory: ${WORK_DIR}"

# Create output directory if it doesn't exist
OUTPUT_DIR="${WORK_DIR}/Result/Filter"
mkdir -p "${OUTPUT_DIR}"

# Create log directory
LOG_DIR="${WORK_DIR}/slurm_logs"
mkdir -p "${LOG_DIR}"

echo "Output will be saved to: ${OUTPUT_DIR}"
echo ""

# Run the analysis
echo "🚀 Starting R Analysis..."
echo "========================"

# Set environment variables to suppress warnings and control threading
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export R_LIBS_USER="${WORK_DIR}/Code/INLA/local_packages"

# Run the main analysis script
echo "Executing BYM_INLA_Production.R..."

Rscript "${WORK_DIR}/Code/INLA/BYM_INLA_Production.R" \
  --config "${WORK_DIR}/Code/INLA/config/analysis_config.yaml" \
  --pesticide-category "compound:${COMPOUND_ID}" \
  --measure-type "${MEASURE_TYPES}" \
  --estimate-types "${ESTIMATE_TYPES}" \
  --lag-years "${LAG_YEARS}" \
  --model-types "${MODEL_TYPES}" \
  --disease-code "${DISEASE_CODE}" \
  --output-file "${OUTPUT_DIR}/Results_${DISEASE_CODE}_${COMPOUND_NAME}_${TIMESTAMP}.csv" \
  --verbose

# Check if the analysis was successful
EXIT_CODE=$?

echo ""
echo "============================================================"

if [ $EXIT_CODE -eq 0 ]; then
    echo "✅ Analysis completed successfully!"
    
    # Check output file
    OUTPUT_FILE="${OUTPUT_DIR}/Results_${DISEASE_CODE}_${COMPOUND_NAME}_${TIMESTAMP}.csv"
    if [ -f "$OUTPUT_FILE" ]; then
        NUM_RESULTS=$(wc -l < "$OUTPUT_FILE")
        echo "📊 Results saved: $((NUM_RESULTS - 1)) records"
        echo "📁 Output file: $OUTPUT_FILE"
        
        # Show a preview of results
        echo ""
        echo "📋 Preview of results (first 5 lines):"
        echo "-------------------------------------"
        head -n 5 "$OUTPUT_FILE"
        
        # Show summary statistics
        echo ""
        echo "📈 Summary Statistics:"
        echo "--------------------"
        echo "Total combinations tested: $((NUM_RESULTS - 1))"
        
        # Count by model type
        echo ""
        echo "Results by Model Type:"
        tail -n +2 "$OUTPUT_FILE" | cut -d',' -f8 | sort | uniq -c | while read count model; do
            echo "  $model: $count results"
        done
        
        # Count by measure type
        echo ""
        echo "Results by Measure Type:"
        tail -n +2 "$OUTPUT_FILE" | cut -d',' -f5 | sort | uniq -c | while read count measure; do
            echo "  $measure: $count results"
        done
        
        # Check for any errors in results
        echo ""
        echo "Error Summary:"
        ERROR_COUNT=$(tail -n +2 "$OUTPUT_FILE" | grep -c "ERROR\|FAILED" || echo "0")
        SUCCESS_COUNT=$(tail -n +2 "$OUTPUT_FILE" | grep -c "Success\|OK" || echo "0")
        echo "  Successful: $SUCCESS_COUNT"
        echo "  Failed: $ERROR_COUNT"
        
    else
        echo "❌ Output file not created: $OUTPUT_FILE"
        EXIT_CODE=1
    fi
else
    echo "❌ Analysis failed with exit code: $EXIT_CODE"
    echo "Check the error messages above for details."
fi

echo ""
echo "Total run time: $((SECONDS / 60)) minutes $((SECONDS % 60)) seconds"
echo "Completed at: $(date)"
echo "============================================================"

exit $EXIT_CODE