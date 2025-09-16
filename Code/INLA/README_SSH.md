# WDP BYM INLA Remote Server Execution Guide

## Overview

This guide documents how to deploy and run the WDP BYM INLA analysis on remote HPC servers via SSH. The system supports automated deployment, execution monitoring, and result retrieval.

## Table of Contents

1. [Server Configuration](#server-configuration)
2. [Initial Setup](#initial-setup)
3. [Deployment Process](#deployment-process)
4. [Running Analyses](#running-analyses)
5. [Monitoring Progress](#monitoring-progress)
6. [Retrieving Results](#retrieving-results)
7. [Troubleshooting](#troubleshooting)
8. [Advanced Usage](#advanced-usage)

---

## Server Configuration

### Connection Details

```bash
Host: cancon.hpccube.com
Port: 65023
Username: acf4pijnzl
Base Path: /public/home/acf4pijnzl
```

### SSH Key Authentication

The server uses SSH key authentication with an expiring RSA key:
- Key file: `acf4pijnzl_cancon.hpccube.com_RsaKeyExpireTime_2025-12-14_15-45-43.txt`
- Expiration: December 14, 2025
- Permissions: Must be set to 600 (`chmod 600 keyfile.txt`)

---

## Initial Setup

### 1. Verify SSH Key

```bash
# Check key exists and has correct permissions
ls -la *.txt | grep RsaKey
chmod 600 acf4pijnzl_cancon.hpccube.com_RsaKeyExpireTime_2025-12-14_15-45-43.txt
```

### 2. Test Connection

```bash
# Test basic SSH connection
ssh -i acf4pijnzl_cancon.hpccube.com_RsaKeyExpireTime_2025-12-14_15-45-43.txt \
    -p 65023 acf4pijnzl@cancon.hpccube.com "echo 'Connection successful'"
```

### 3. Set Up SSH Alias (Optional)

Add to `~/.ssh/config`:

```
Host wdp-hpc
    HostName cancon.hpccube.com
    Port 65023
    User acf4pijnzl
    IdentityFile ~/path/to/acf4pijnzl_cancon.hpccube.com_RsaKeyExpireTime_2025-12-14_15-45-43.txt
    StrictHostKeyChecking no
```

Then connect simply with: `ssh wdp-hpc`

---

## Deployment Process

### Automated Deployment Script

The primary deployment script handles all aspects of remote setup:

```bash
./Code/INLA/deploy_to_server_fixed.sh
```

This script performs the following operations:

1. **Creates Remote Directory Structure**
   ```
   /public/home/acf4pijnzl/
   ├── WDP_Analysis/
   │   ├── Code/INLA/
   │   ├── Data/Processed/
   │   └── Results/
   └── WDP_Results/
       └── {analysis_name}/
   ```

2. **Transfers Code Files**
   - Main analysis script (`BYM_INLA_Production.R`)
   - Utility modules (`utils/*.R`)
   - Configuration files (`config/*.yaml`)

3. **Handles Data Transfer**
   - Checks if data already exists on server
   - Transfers only if necessary (~5GB compressed)
   - Preserves directory structure

4. **Creates Execution Scripts**
   - Analysis runner with proper paths
   - Configuration modifications
   - Environment setup

### Manual Deployment (Advanced)

For custom deployments:

```bash
# 1. Create remote directories
ssh -i $SSH_KEY -p 65023 acf4pijnzl@cancon.hpccube.com \
    "mkdir -p /public/home/acf4pijnzl/WDP_Analysis/{Code,Data,Results}"

# 2. Transfer code
rsync -avz -e "ssh -i $SSH_KEY -p 65023" \
    Code/INLA/ acf4pijnzl@cancon.hpccube.com:/public/home/acf4pijnzl/WDP_Analysis/Code/INLA/

# 3. Transfer data (if needed)
tar -czf - -C Data/Processed . | \
    ssh -i $SSH_KEY -p 65023 acf4pijnzl@cancon.hpccube.com \
    "cd /public/home/acf4pijnzl/WDP_Analysis/Data/Processed && tar -xzf -"
```

---

## Running Analyses

### Standard Analysis Workflow

1. **Deploy and Start Analysis**
   ```bash
   # Run deployment script
   ./Code/INLA/deploy_to_server_fixed.sh
   
   # Note the analysis name (e.g., 24D_Full_Analysis_20240915_143022)
   ```

2. **Analysis Parameters**
   
   The deployment script creates a server-side runner that executes:
   - Multiple exposure measures (Weight, Density)
   - Multiple estimates (min, avg, max)
   - Multiple models (M0, M1, M2, M3)
   - Multiple dose-response types (Linear, Non-linear)

### Direct Server Execution

For custom analyses, SSH into the server and run directly:

```bash
# Connect to server
ssh -i $SSH_KEY -p 65023 acf4pijnzl@cancon.hpccube.com

# Navigate to project
cd /public/home/acf4pijnzl/WDP_Analysis

# Run analysis
Rscript Code/INLA/BYM_INLA_Production.R \
    --pesticide-category compound:2 \
    --measure-type Weight,Density \
    --estimate-types min,avg,max \
    --lag-years 5 \
    --model-types M0,M1,M2,M3 \
    --disease-code C81-C96 \
    --output-file Results/my_analysis.csv \
    --verbose
```

### Using nohup for Long-Running Analyses

```bash
nohup Rscript Code/INLA/BYM_INLA_Production.R \
    --pesticide-category ALL \
    --measure-type Weight,Density \
    --estimate-types min,avg,max \
    --lag-years 5,10 \
    --model-types M0,M1,M2,M3 \
    --disease-code C81-C96 \
    > analysis.log 2>&1 &

# Save PID for monitoring
echo $! > analysis.pid
```

---

## Monitoring Progress

### Using the Monitor Script

After deployment, use the generated monitoring script:

```bash
./monitor_analysis.sh
```

This script shows:
- Process status (running/completed)
- Recent log output (last 30 lines)
- Result files created
- Current progress

### Manual Monitoring

```bash
# Check if analysis is running
ssh -i $SSH_KEY -p 65023 acf4pijnzl@cancon.hpccube.com \
    "ps aux | grep BYM_INLA | grep -v grep"

# View recent log output
ssh -i $SSH_KEY -p 65023 acf4pijnzl@cancon.hpccube.com \
    "tail -n 50 /public/home/acf4pijnzl/WDP_Results/*/server_execution.log"

# Check result files
ssh -i $SSH_KEY -p 65023 acf4pijnzl@cancon.hpccube.com \
    "ls -la /public/home/acf4pijnzl/WDP_Results/*/Results_*.csv"
```

### Real-time Log Streaming

```bash
# Stream logs in real-time
ssh -i $SSH_KEY -p 65023 acf4pijnzl@cancon.hpccube.com \
    "tail -f /public/home/acf4pijnzl/WDP_Results/*/server_execution.log"
```

---

## Retrieving Results

### Using the Download Script

After analysis completion:

```bash
./download_results.sh
```

This downloads all results to: `Results/{analysis_name}/`

### Manual Download

```bash
# Download specific files
scp -i $SSH_KEY -P 65023 \
    acf4pijnzl@cancon.hpccube.com:/public/home/acf4pijnzl/WDP_Results/*/Results_*.csv \
    Results/

# Download entire analysis directory
rsync -avz -e "ssh -i $SSH_KEY -p 65023" \
    acf4pijnzl@cancon.hpccube.com:/public/home/acf4pijnzl/WDP_Results/{analysis_name}/ \
    Results/{analysis_name}/
```

### Result Files

Typical analysis produces:
- `Results_24D_Linear.csv` - Linear dose-response models
- `Results_24D_Mixed.csv` - Mixed linear/non-linear models
- `Results_24D_Nonlinear.csv` - Non-linear dose-response models
- `Analysis_Summary.txt` - Summary statistics
- `*.log` - Execution logs

---

## Troubleshooting

### Common Issues and Solutions

#### 1. SSH Connection Failures

```bash
# Issue: Permission denied
Solution: chmod 600 on SSH key file

# Issue: Connection timeout
Solution: Check firewall, verify port 65023 is open

# Issue: Host key verification failed
Solution: Add StrictHostKeyChecking no to SSH command
```

#### 2. Path Issues

```bash
# Issue: "Directory not found"
# Server uses /public/home not /home
Correct: /public/home/acf4pijnzl/
Wrong: /home/acf4pijnzl/
```

#### 3. R Package Issues

```bash
# Connect and install missing packages
ssh -i $SSH_KEY -p 65023 acf4pijnzl@cancon.hpccube.com

R
install.packages("INLA", repos=c(getOption("repos"),
                 INLA="https://inla.r-inla-download.org/R/stable"))
q()
```

#### 4. Memory/Resource Issues

```bash
# Check server resources
ssh -i $SSH_KEY -p 65023 acf4pijnzl@cancon.hpccube.com \
    "free -h; df -h /public/home/acf4pijnzl"

# Monitor during execution
ssh -i $SSH_KEY -p 65023 acf4pijnzl@cancon.hpccube.com "htop"
```

### Debug Mode Execution

For troubleshooting, run with verbose output:

```bash
ssh -i $SSH_KEY -p 65023 acf4pijnzl@cancon.hpccube.com
cd /public/home/acf4pijnzl/WDP_Analysis

# Run single model for testing
Rscript Code/INLA/BYM_INLA_Production.R \
    --pesticide-category TEST \
    --measure-type Weight \
    --estimate-types avg \
    --lag-years 5 \
    --model-types M0 \
    --disease-code C81-C96 \
    --verbose \
    --dry-run  # Validates without running
```

---

## Advanced Usage

### Parallel Execution

Run multiple analyses simultaneously:

```bash
# Create batch script on server
cat > run_parallel.sh << 'EOF'
#!/bin/bash
# Run multiple compounds in parallel

compounds=(2 5 10 15 20)

for comp in "${compounds[@]}"; do
    nohup Rscript Code/INLA/BYM_INLA_Production.R \
        --pesticide-category compound:$comp \
        --measure-type Weight,Density \
        --estimate-types min,avg,max \
        --lag-years 5 \
        --model-types M0,M1,M2,M3 \
        --disease-code C81-C96 \
        > compound_${comp}.log 2>&1 &
    
    echo "Started compound $comp with PID $!"
    sleep 5  # Stagger starts
done
EOF

chmod +x run_parallel.sh
./run_parallel.sh
```

### Custom Configuration

Modify analysis configuration on server:

```bash
# Edit configuration remotely
ssh -i $SSH_KEY -p 65023 acf4pijnzl@cancon.hpccube.com \
    "nano /public/home/acf4pijnzl/WDP_Analysis/Code/INLA/config/analysis_config.yaml"

# Or upload local config
scp -i $SSH_KEY -P 65023 \
    my_custom_config.yaml \
    acf4pijnzl@cancon.hpccube.com:/public/home/acf4pijnzl/WDP_Analysis/Code/INLA/config/
```

### Server Maintenance

```bash
# Clean up old results
ssh -i $SSH_KEY -p 65023 acf4pijnzl@cancon.hpccube.com \
    "find /public/home/acf4pijnzl/WDP_Results -name '*.log' -mtime +30 -delete"

# Archive completed analyses
ssh -i $SSH_KEY -p 65023 acf4pijnzl@cancon.hpccube.com \
    "cd /public/home/acf4pijnzl/WDP_Results && \
     tar -czf archive_$(date +%Y%m%d).tar.gz 24D_Full_Analysis_*"
```

---

## Security Best Practices

1. **SSH Key Management**
   - Store key in secure location
   - Set restrictive permissions (600)
   - Note expiration date (2025-12-14)
   - Request new key before expiration

2. **Data Security**
   - Use encrypted connections only
   - Don't store sensitive data in logs
   - Clean up temporary files after analysis

3. **Access Control**
   - Use specific paths, not wildcards
   - Verify file permissions on server
   - Log all analysis runs for audit trail

---

## Performance Tips

1. **Data Transfer Optimization**
   - Use compression for large transfers
   - Transfer only changed files with rsync
   - Keep frequently used data on server

2. **Analysis Optimization**
   - Run during off-peak hours if possible
   - Use appropriate model complexity
   - Monitor resource usage

3. **Result Management**
   - Download results promptly
   - Clean up server space regularly
   - Use compression for archives

---

## Quick Reference Card

```bash
# Deploy and run
./Code/INLA/deploy_to_server_fixed.sh

# Monitor progress
./monitor_analysis.sh

# Download results  
./download_results.sh

# Direct SSH
ssh -i acf4pijnzl_*.txt -p 65023 acf4pijnzl@cancon.hpccube.com

# Key paths
Project: /public/home/acf4pijnzl/WDP_Analysis
Results: /public/home/acf4pijnzl/WDP_Results
```

---

## Support

For issues specific to the HPC server:
- Check server status page (if available)
- Contact HPC support with job ID
- Review `/public/home/acf4pijnzl/WDP_Results/*/server_execution.log`

For WDP analysis issues:
- See main documentation: `Code/INLA/BYM.md`
- Check test instructions: `Code/INLA/docs/Testing_Instructions.md`
- Review code modifications: `Code/INLA/docs/WONDER_Code_Modifications_Summary.md`

**Last Updated**: 2024-09-15  
**Server Key Expiration**: 2025-12-14