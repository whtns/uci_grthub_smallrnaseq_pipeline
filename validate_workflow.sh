#!/bin/bash

# Quick validation script for the Snakemake workflow
# This script performs basic checks before running the full workflow

echo "=== Small RNA Workflow Validation ==="
echo "Date: $(date)"
echo

# Check if required files exist
echo "1. Checking required files..."

if [ -f "Snakefile" ]; then
    echo "   ✓ Snakefile found"
else
    echo "   ✗ Snakefile not found"
    exit 1
fi

if [ -f "config.yaml" ]; then
    echo "   ✓ config.yaml found"
else
    echo "   ✗ config.yaml not found"
    exit 1
fi

if [ -f "cluster.yaml" ]; then
    echo "   ✓ cluster.yaml found"
else
    echo "   ✗ cluster.yaml not found"
    exit 1
fi

# Check if input directory exists
FASTQ_DIR="FastqFiles/hts.igb.uci.edu/rob25092453"
if [ -d "$FASTQ_DIR" ]; then
    echo "   ✓ FASTQ directory found: $FASTQ_DIR"
    
    # Count FASTQ files
    FASTQ_COUNT=$(find "$FASTQ_DIR" -name "*Sequences.txt.gz" | wc -l)
    echo "   ✓ Found $FASTQ_COUNT FASTQ files"
    
    if [ $FASTQ_COUNT -eq 0 ]; then
        echo "   ⚠  Warning: No FASTQ files found matching pattern '*Sequences.txt.gz'"
    fi
else
    echo "   ✗ FASTQ directory not found: $FASTQ_DIR"
    exit 1
fi

echo

# Check for Snakemake
echo "2. Checking software requirements..."
if command -v snakemake &> /dev/null; then
    echo "   ✓ Snakemake available: $(snakemake --version)"
else
    echo "   ✗ Snakemake not found. Try: module load snakemake"
fi

# Test workflow syntax
echo

echo "3. Testing workflow syntax..."
if snakemake --dry-run --quiet > /dev/null 2>&1; then
    echo "   ✓ Workflow syntax is valid"
else
    echo "   ✗ Workflow syntax error. Run 'snakemake --dry-run' for details"
    exit 1
fi

echo

# Show sample information
echo "4. Sample information:"
echo "   Discovering samples from FASTQ files..."

python3 << 'EOF'
import glob
import os

fastq_files = glob.glob("FastqFiles/hts.igb.uci.edu/rob25092453/*Sequences.txt.gz")

samples = {}
for fastq_file in fastq_files:
    filename = os.path.basename(fastq_file)
    parts = filename.replace("-Sequences.txt.gz", "").split("-")
    read_part = parts[-1]
    sample_base = "-".join(parts[:-1])
    
    if sample_base not in samples:
        samples[sample_base] = {}
    samples[sample_base][read_part] = fastq_file

paired_samples = {k: v for k, v in samples.items() if "READ1" in v and "READ2" in v}

print(f"   Found {len(paired_samples)} paired samples:")
for i, sample in enumerate(sorted(paired_samples.keys()), 1):
    print(f"   {i:2d}. {sample}")

if len(paired_samples) == 0:
    print("   ⚠  Warning: No paired samples found")
EOF

echo

echo "5. Recommended next steps:"
echo "   - For dry run: snakemake --dry-run"
echo "   - For local execution: snakemake --cores 4"
echo "   - For cluster execution: sbatch submit_workflow.sh"
echo

echo "=== Validation Complete ==="