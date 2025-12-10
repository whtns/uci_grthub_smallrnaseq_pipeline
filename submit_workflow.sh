#!/bin/bash
#SBATCH --job-name=snakemake_smallRNA
#SBATCH -A SBSANDME_LAB
#SBATCH -p standard
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=08:00:00
#SBATCH --error=snakemake_master.%j.err
#SBATCH --output=snakemake_master.%j.out

# Small RNA processing workflow submission script
# This script submits the Snakemake workflow to SLURM

# Load required modules
module load python/3.8
module load snakemake

# Create necessary directories
mkdir -p logs

# Print workflow information
echo "Starting Small RNA processing workflow"
echo "Date: $(date)"
echo "Working directory: $(pwd)"
echo "Snakemake version: $(snakemake --version)"

# Dry run first to check workflow
echo "Performing dry run..."
snakemake --dry-run --quiet

if [ $? -eq 0 ]; then
    echo "Dry run successful. Starting workflow execution..."
    
    # Run the workflow with SLURM cluster support
    snakemake \
        --cluster-config cluster.yaml \
        --cluster "sbatch --account={cluster.account} --partition={cluster.partition} --nodes={cluster.nodes} --ntasks={cluster.ntasks} --cpus-per-task={cluster.cpus} --mem={cluster.mem} --time={cluster.time} --output={cluster.output} --error={cluster.error}" \
        --jobs 10 \
        --latency-wait 30 \
        --keep-going \
        --rerun-incomplete
        
    echo "Workflow completed with exit code: $?"
    echo "End time: $(date)"
    
    # Generate summary
    snakemake create_summary
    
else
    echo "Dry run failed. Please check the workflow configuration."
    exit 1
fi