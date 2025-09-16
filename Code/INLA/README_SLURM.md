# WDP SLURM Job Submission Guide

## Overview

This guide explains how to run WDP analyses on HPC clusters using the SLURM workload manager. Two job scripts are provided:

1. `HCP_RUN_Test1.sh` - A test job that runs a single model to verify the setup
2. `HCP_RUN_Full_24D.sh` - A full 2,4-D analysis with all model combinations

## Prerequisites

1. You must be on the HPC server
2. The WDP project must be deployed to `/public/home/acf4pijnzl/WDP_Analysis`
3. Data files must be available in `/public/home/acf4pijnzl/WDP_Analysis/Data/Processed`

## SLURM Job Scripts

### HCP_RUN_Test1.sh

This script runs a minimal test analysis to verify that everything is working correctly:

- Runs a single model (M0) with avg estimate for Weight measure
- Uses the TEST pesticide category
- Completes relatively quickly (typically within 30 minutes)

### HCP_RUN_Full_24D.sh

This script runs the complete 2,4-D analysis:

- All measure types: Weight, Density
- All estimate types: min, avg, max
- All model types: M0, M1, M2, M3
- Both linear and non-linear dose-response models
- Expected to take 24-48 hours to complete

## Resource Requirements

Both scripts are configured with the following resources:

```
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=64G
#SBATCH --time=48:00:00
```

## Job Submission

### Method 1: Using the submission script

```bash
# Submit the test job
./Code/INLA/submit_slurm_job.sh test

# Submit the full analysis job
./Code/INLA/submit_slurm_job.sh full
```

### Method 2: Direct submission

```bash
# Submit the test job
sbatch Code/INLA/HCP_RUN_Test1.sh

# Submit the full analysis job
sbatch Code/INLA/HCP_RUN_Full_24D.sh
```

## Monitoring Jobs

### Check job status

```bash
# View all your jobs
squeue -u acf4pijnzl

# View specific job
squeue -j <job_id>
```

### View job output

```bash
# View output file (updated in real-time)
tail -f /public/home/acf4pijnzl/WDP_Results/slurm-<job_id>.out

# View error file
tail -f /public/home/acf4pijnzl/WDP_Results/slurm-<job_id>.err
```

## Results

Results are saved to timestamped directories in:
```
/public/home/acf4pijnzl/WDP_Results/
```

Each job creates a directory with the format:
- Test jobs: `Test_Analysis_YYYYMMDD_HHMMSS`
- Full analysis: `24D_Full_Analysis_YYYYMMDD_HHMMSS`

## Customization

To modify the resource requirements, edit the SLURM directives at the top of each script:

```bash
#SBATCH --partition=compute     # Partition name
#SBATCH --nodes=1               # Number of nodes
#SBATCH --ntasks-per-node=8     # CPU cores per node
#SBATCH --mem=64G               # Memory allocation
#SBATCH --time=48:00:00         # Time limit (HH:MM:SS)
```

To modify analysis parameters, edit the variables in the script:

```bash
COMPOUND_ID="2"
COMPOUND_NAME="2,4-D"
DISEASE_CODE="C81-C96"
MEASURE_TYPES="Weight,Density"
ESTIMATE_TYPES="min,avg,max"
LAG_YEARS="5"
MODEL_TYPES="M0,M1,M2,M3"
```

## Troubleshooting

### Common Issues

1. **Module loading fails**
   - The script attempts to load R module automatically
   - If it fails, you may need to load it manually:
     ```bash
     module load R/4.3.0
     ```

2. **Package installation issues**
   - The script checks and installs required packages automatically
   - If installation fails, try running in an interactive session:
     ```bash
     R
     install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"))
     ```

3. **Permission errors**
   - Ensure the scripts have execute permissions:
     ```bash
     chmod +x Code/INLA/HCP_RUN_Test1.sh
     chmod +x Code/INLA/HCP_RUN_Full_24D.sh
     chmod +x Code/INLA/submit_slurm_job.sh
     ```

### Canceling Jobs

```bash
# Cancel a specific job
scancel <job_id>

# Cancel all your jobs
scancel -u acf4pijnzl
```

## Support

For issues with the SLURM scripts, contact the WDP development team.

For HPC-specific issues, contact your system administrator.