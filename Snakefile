#!/usr/bin/env python3

"""
Snakemake workflow for processing small RNA sequencing data
Based on the existing shell scripts in shellscripts/

This workflow:
1. Runs FastQC on raw FASTQ files
2. Trims adapters using Cutadapt for small RNA v4 kit
3. Runs FastQC on trimmed files (optional)
"""

import os
import glob
from pathlib import Path

# Configuration
configfile: "config.yaml"

# Define base paths
FASTQ_DIR = "data/FASTQ"
FASTQC_DIR = "fastQC"
TRIMMED_DIR = "trimmed_fastq"
LOGS_DIR = "logs"

# Get all FASTQ files matching the pattern
FASTQ_FILES = glob.glob(f"{FASTQ_DIR}/*Sequences.txt.gz")

# Extract sample information from filenames
def get_sample_info(fastq_path):
    """Extract sample name and read info from FASTQ filename"""
    filename = os.path.basename(fastq_path)
    # Example: xR064-L2-G4-P18-AAGGCCTG-TCTCAATT-READ1-Sequences.txt.gz
    parts = filename.replace("-Sequences.txt.gz", "").split("-")
    read_part = parts[-1]  # READ1 or READ2
    sample_base = "-".join(parts[:-1])  # Everything except READx
    return sample_base, read_part

# Get unique samples (paired-end)
SAMPLES = {}
for fastq_file in FASTQ_FILES:
    sample_base, read = get_sample_info(fastq_file)
    if sample_base not in SAMPLES:
        SAMPLES[sample_base] = {}
    SAMPLES[sample_base][read] = fastq_file

# Filter to only include samples with both READ1 and READ2
PAIRED_SAMPLES = {k: v for k, v in SAMPLES.items() 
                  if "READ1" in v and "READ2" in v}

print(f"Found {len(PAIRED_SAMPLES)} paired samples:")
for sample in PAIRED_SAMPLES:
    print(f"  - {sample}")

# Define all output files
rule all:
    input:
        # FastQC outputs for raw files
        expand("{fastqc_dir}/html/{sample}-{read}-Sequences_fastqc.html",
               fastqc_dir=FASTQC_DIR,
               sample=PAIRED_SAMPLES.keys(),
               read=["READ1", "READ2"]),
        # Cutadapt trimmed outputs
        expand("{trimmed_dir}/{sample}.cut.R{read_num}.fastq",
               trimmed_dir=TRIMMED_DIR,
               sample=PAIRED_SAMPLES.keys(),
               read_num=[1, 2]),
        # # FastQC outputs for trimmed files (optional)
        # expand("{fastqc_dir}/trimmed_html/{sample}.cut.R{read_num}_fastqc.html",
        #        fastqc_dir=FASTQC_DIR,
        #        sample=PAIRED_SAMPLES.keys(),
        #        read_num=[1, 2])

# Rule to run FastQC on raw FASTQ files
rule fastqc_raw:
    input:
        r1 = lambda wildcards: PAIRED_SAMPLES[wildcards.sample]["READ1"],
        r2 = lambda wildcards: PAIRED_SAMPLES[wildcards.sample]["READ2"]
    output:
        html_r1 = FASTQC_DIR + "/html/{sample}-READ1-Sequences_fastqc.html",
        html_r2 = FASTQC_DIR + "/html/{sample}-READ2-Sequences_fastqc.html",
        zip_r1 = FASTQC_DIR + "/{sample}-READ1-Sequences_fastqc.zip",
        zip_r2 = FASTQC_DIR + "/{sample}-READ2-Sequences_fastqc.zip"
    params:
        outdir = FASTQC_DIR,
        html_dir = FASTQC_DIR + "/html"
    threads: 8
    resources:
        mem_mb = 4000,
        runtime = 60
    log:
        LOGS_DIR + "/fastqc_raw/{sample}.log"
    shell:
        """
        module load fastqc/0.11.9
        
        # Create output directories
        mkdir -p {params.outdir}
        mkdir -p {params.html_dir}
        
        # Run FastQC
        fastqc -t {threads} {input.r1} {input.r2} -o {params.outdir} 2>&1 | tee {log}
        
        # Move HTML files to html subdirectory
        mv {params.outdir}/*-READ*-Sequences_fastqc.html {params.html_dir}/
        
        module unload fastqc/0.11.9
        """

# Rule to trim adapters using Cutadapt
rule cutadapt_trim:
    input:
        r1 = lambda wildcards: PAIRED_SAMPLES[wildcards.sample]["READ1"],
        r2 = lambda wildcards: PAIRED_SAMPLES[wildcards.sample]["READ2"]
    output:
        trimmed_r1 = TRIMMED_DIR + "/{sample}.cut.R1.fastq.gz",
        trimmed_r2 = TRIMMED_DIR + "/{sample}.cut.R2.fastq.gz"
    params:
        adapter_r1 = "TGGAATTCTCGGGTGCCAAGG",
        adapter_r2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGA",
        min_length = 16,
        quality_cutoff = 20,
        trimmed_dir = TRIMMED_DIR
    threads: 4
    resources:
        mem_mb = 8000,
        runtime = 120
    log:
        LOGS_DIR + "/cutadapt/{sample}.log"
    shell:
        """
        module load cutadapt/2.10
        
        # Create output directory
        mkdir -p {params.trimmed_dir}
        
        # Trim adapters using cutadapt
        # Biooscientific small RNA v4 kit requires special trimming
        cutadapt --pair-adapters \\
            --quality-cutoff {params.quality_cutoff} \\
            --adapter {params.adapter_r1} \\
            -A {params.adapter_r2} \\
            --output {output.trimmed_r1} \\
            --paired-output {output.trimmed_r2} \\
            --minimum-length {params.min_length} \\
            --cores {threads} \\
            {input.r1} {input.r2} > {log} 2>&1
        
        module unload cutadapt/2.10
        """

# Rule to run FastQC on trimmed files
rule fastqc_trimmed:
    input:
        trimmed_r1 = TRIMMED_DIR + "/{sample}.cut.R1.fastq.gz",
        trimmed_r2 = TRIMMED_DIR + "/{sample}.cut.R2.fastq.gz"
    output:
        html_r1 = FASTQC_DIR + "/trimmed_html/{sample}.cut.R1_fastqc.html",
        html_r2 = FASTQC_DIR + "/trimmed_html/{sample}.cut.R2_fastqc.html",
        zip_r1 = FASTQC_DIR + "/trimmed/{sample}.cut.R1_fastqc.zip",
        zip_r2 = FASTQC_DIR + "/trimmed/{sample}.cut.R2_fastqc.zip"
    params:
        outdir = FASTQC_DIR + "/trimmed",
        html_dir = FASTQC_DIR + "/trimmed_html"
    threads: 8
    resources:
        mem_mb = 4000,
        runtime = 60
    log:
        LOGS_DIR + "/fastqc_trimmed/{sample}.log"
    shell:
        """
        module load fastqc/0.11.9
        
        # Create output directories
        mkdir -p {params.outdir}
        mkdir -p {params.html_dir}
        
        # Run FastQC on trimmed files
        fastqc -t {threads} {input.trimmed_r1} {input.trimmed_r2} -o {params.outdir} 2>&1 | tee {log}
        
        # Move HTML files to html subdirectory
        mv {params.outdir}/*.html {params.html_dir}/
        
        module unload fastqc/0.11.9
        """

# Rule to create a summary report
rule create_summary:
    input:
        expand("{trimmed_dir}/{sample}.cut.R{read_num}.fastq",
               trimmed_dir=TRIMMED_DIR,
               sample=PAIRED_SAMPLES.keys(),
               read_num=[1, 2])
    output:
        "processing_summary.txt"
    params:
        sample_list = list(PAIRED_SAMPLES.keys()),
        fastqc_dir = FASTQC_DIR,
        trimmed_dir = TRIMMED_DIR,
        logs_dir = LOGS_DIR
    shell:
        """
        echo "Small RNA Processing Summary" > {output}
        echo "Generated on: $(date)" >> {output}
        echo "=========================" >> {output}
        echo "" >> {output}
        echo "Processed samples:" >> {output}
        for sample in {params.sample_list}; do
            echo "  - $sample" >> {output}
        done
        echo "" >> {output}
        echo "Output files:" >> {output}
        echo "  - Raw FastQC reports: {params.fastqc_dir}/html/" >> {output}
        echo "  - Trimmed FASTQ files: {params.trimmed_dir}/" >> {output}
        echo "  - Trimmed FastQC reports: {params.fastqc_dir}/trimmed_html/" >> {output}
        echo "  - Log files: {params.logs_dir}/" >> {output}
        """

# Utility rules for specific processing steps
rule fastqc_only:
    input:
        expand("{fastqc_dir}/html/{sample}-{read}-Sequences_fastqc.html",
               fastqc_dir=FASTQC_DIR,
               sample=PAIRED_SAMPLES.keys(),
               read=["READ1", "READ2"])

rule cutadapt_only:
    input:
        expand("{trimmed_dir}/{sample}.cut.R{read_num}.fastq",
               trimmed_dir=TRIMMED_DIR,
               sample=PAIRED_SAMPLES.keys(),
               read_num=[1, 2])