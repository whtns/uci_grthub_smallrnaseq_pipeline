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
FASTQ_FILES = glob.glob(f"{FASTQ_DIR}/*fastq.gz")

# Extract sample information from filenames
def get_sample_info(fastq_path):
    """Extract sample name and read info from FASTQ filename"""
    filename = os.path.basename(fastq_path)
    # Example: xR064-L2-G4-P18-AAGGCCTG-TCTCAATT-R1.fastq.gz
    parts = filename.replace(".fastq.gz", "").split("-")
    read_part = parts[-1]  # R1 or R2
    sample_base = "-".join(parts[:-1])  # Everything except Rx
    return sample_base, read_part

# Get unique samples - R1 only
SAMPLES = {}
for fastq_file in FASTQ_FILES:
    sample_base, read = get_sample_info(fastq_file)
    if sample_base not in SAMPLES:
        SAMPLES[sample_base] = {}
    SAMPLES[sample_base][read] = fastq_file

# Filter to only include samples with R1
PAIRED_SAMPLES = {k: v["R1"] for k, v in SAMPLES.items() 
                  if "R1" in v}

print(f"Found {len(PAIRED_SAMPLES)} R1 samples:")
for sample in PAIRED_SAMPLES:
    print(f"  - {sample}")

# Define all output files
rule all:
    input:
        # FastQC outputs for raw R1 files
        expand("{fastqc_dir}/html/{sample}-R1_fastqc.html",
               fastqc_dir=FASTQC_DIR,
               sample=PAIRED_SAMPLES.keys()),
        # Cutadapt trimmed R1 outputs
        expand("{trimmed_dir}/{sample}.cut.R1.fastq",
               trimmed_dir=TRIMMED_DIR,
               sample=PAIRED_SAMPLES.keys()),
        # # FastQC outputs for trimmed files (optional)
        # expand("{fastqc_dir}/trimmed_html/{sample}.cut.R{read_num}_fastqc.html",
        #        fastqc_dir=FASTQC_DIR,
        #        sample=PAIRED_SAMPLES.keys(),
        #        read_num=[1, 2])

# Rule to run FastQC on raw R1 FASTQ files
rule fastqc_raw:
    input:
        r1 = lambda wildcards: PAIRED_SAMPLES[wildcards.sample]
    output:
        html_r1 = FASTQC_DIR + "/html/{sample}-R1_fastqc.html",
        zip_r1 = FASTQC_DIR + "/{sample}-R1_fastqc.zip"
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
        
        # Run FastQC on R1 only
        fastqc -t {threads} {input.r1} -o {params.outdir} 2>&1 | tee {log}
        
        # Move HTML files to html subdirectory
        mv {params.outdir}/*-R1_fastqc.html {params.html_dir}/
        
        module unload fastqc/0.11.9
        """

# Rule to trim adapters using Cutadapt on R1 only
rule cutadapt_trim:
    input:
        r1 = lambda wildcards: PAIRED_SAMPLES[wildcards.sample]
    output:
        trimmed_r1 = TRIMMED_DIR + "/{sample}.cut.R1.fastq"
    params:
        adapter_1 = "CTGTCTCTTATACACATCT",
        adapter_2 = "TGGAATTCTCGGGTGCCAAGG",
        min_length = 15,
        trim_5 = 4,
        trim_3 = 4,
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
        cutadapt -a {params.adapter_1} \
            -a {params.adapter_2} \
            -m {params.min_length} \
            -u {params.trim_5} \
            -u -{params.trim_3} \
            --output {output.trimmed_r1} \
            --cores {threads} \
            {input.r1} > {log} 2>&1
        
        module unload cutadapt/2.10
        """

# Rule to run FastQC on trimmed R1 files
rule fastqc_trimmed:
    input:
        trimmed_r1 = TRIMMED_DIR + "/{sample}.cut.R1.fastq"
    output:
        html_r1 = FASTQC_DIR + "/trimmed_html/{sample}.cut.R1_fastqc.html",
        zip_r1 = FASTQC_DIR + "/trimmed/{sample}.cut.R1_fastqc.zip"
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
        
        # Run FastQC on trimmed R1 files
        fastqc -t {threads} {input.trimmed_r1} -o {params.outdir} 2>&1 | tee {log}
        
        # Move HTML files to html subdirectory
        mv {params.outdir}/*.html {params.html_dir}/
        
        module unload fastqc/0.11.9
        """

# Rule to create a summary report
rule create_summary:
    input:
        expand("{trimmed_dir}/{sample}.cut.R1.fastq",
               trimmed_dir=TRIMMED_DIR,
               sample=PAIRED_SAMPLES.keys())
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
        expand("{fastqc_dir}/html/{sample}-R1_fastqc.html",
               fastqc_dir=FASTQC_DIR,
               sample=PAIRED_SAMPLES.keys())

rule cutadapt_only:
    input:
        expand("{trimmed_dir}/{sample}.cut.R1.fastq",
               trimmed_dir=TRIMMED_DIR,
               sample=PAIRED_SAMPLES.keys())