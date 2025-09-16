#!/bin/bash

# Script to submit SLURM jobs for WDP analysis
# Usage: ./submit_slurm_job.sh [test|full]

echo "WDP SLURM Job Submission Script"
echo "==============================="

# Check if we're on the HPC server
if [[ ! -d "/public/home" ]]; then
    echo "❌ Error: This script must be run on the HPC server"
    echo "Please SSH to the server first:"
    echo "ssh -i acf4pijnzl_cancon.hpccube.com_RsaKeyExpireTime_2025-12-14_15-45-43.txt -p 65023 acf4pijnzl@cancon.hpccube.com"
    exit 1
fi

# Determine which job to submit
JOB_TYPE="full"
if [[ $# -gt 0 ]]; then
    JOB_TYPE="$1"
fi

# Set working directory
cd /public/home/acf4pijnzl/WDP_Analysis

echo "Current directory: $(pwd)"

# Submit the appropriate job
if [[ "$JOB_TYPE" == "test" ]]; then
    echo "Submitting test job..."
    sbatch Code/INLA/HCP_RUN_Test1.sh
elif [[ "$JOB_TYPE" == "full" ]]; then
    echo "Submitting full 2,4-D analysis job..."
    sbatch Code/INLA/HCP_RUN_Full_24D.sh
else
    echo "❌ Unknown job type: $JOB_TYPE"
    echo "Usage: ./submit_slurm_job.sh [test|full]"
    exit 1
fi

if [[ $? -eq 0 ]]; then
    echo "✅ Job submitted successfully!"
    echo "Use 'squeue -u acf4pijnzl' to monitor job status"
else
    echo "❌ Failed to submit job"
    exit 1
fi