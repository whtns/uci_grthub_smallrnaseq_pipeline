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
# Select the library-prep kit by pointing at the matching sub-config:
#   config_v3.yaml  -> Bioo Scientific NEXTFLEX Small RNA-Seq v3 (4N randomers)
#   config_v4.yaml  -> Revvity NEXTFLEX Small RNA-Seq v4 (no end trimming)
configfile: "config_v3.yaml"

# Define base paths
# Inputs live under data/; all generated outputs go under OUTPUT_DIR.
FASTQ_DIR = "data/FASTQ"
OUTPUT_DIR = "output"
FASTQC_DIR = OUTPUT_DIR + "/fastQC"
TRIMMED_DIR = OUTPUT_DIR + "/trimmed_fastq"
LOGS_DIR = OUTPUT_DIR + "/logs"

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

print(f"Found {len(PAIRED_SAMPLES)} R1 samples:", file=sys.stderr)
for sample in PAIRED_SAMPLES:
    print(f"  - {sample}", file=sys.stderr)

# Define all output files
rule all:
    input:
        OUTPUT_DIR + "/multiqc_report.html",
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

# Rule to trim adapters using Cutadapt on R1 only (v3 - commented out, use v4 instead)
# rule cutadapt_trim_v3:
#     input:
#         r1 = lambda wildcards: PAIRED_SAMPLES[wildcards.sample]
#     output:
#         trimmed_r1 = TRIMMED_DIR + "/{sample}.cut.R1.fastq"
#     params:
#         adapter_1 = "CTGTCTCTTATACACATCT",
#         adapter_2 = "TGGAATTCTCGGGTGCCAAGG",
#         min_length = 15,
#         trim_5 = 4,
#         trim_3 = 4,
#         trimmed_dir = TRIMMED_DIR
#     threads: 4
#     resources:
#         mem_mb = 24000
#     log:
#         LOGS_DIR + "/cutadapt/{sample}.log"
#     shell:
#         """
#         module load cutadapt/2.10
#         
#         # Create output directory
#         mkdir -p {params.trimmed_dir}
#         
#         # Trim adapters using cutadapt
#         cutadapt -a {params.adapter_1} \
#             -a {params.adapter_2} \
#             -m {params.min_length} \
#             -u {params.trim_5} \
#             -u -{params.trim_3} \
#             --output {output.trimmed_r1} \
#             --cores {threads} \
#             {input.r1} > {log} 2>&1
#         
#         module unload cutadapt/2.10
#         """

# Rule to trim adapters using Cutadapt on R1 only.
# Kit-specific behavior is driven by the active configfile (config_v3 / config_v4):
#   - adapter, min_length, quality_cutoff
#   - trim_5 / trim_3: fixed bases removed from each end AFTER adapter trimming
#     (v3 = 4/4 for the 4N randomized adapters; v4 = 0/0)
# IMPORTANT: cutadapt applies -u (--cut) BEFORE adapter trimming, so for the v3
# randomers we must trim the adapter first, then remove the flanking bases in a
# second pass; otherwise -u -N would cut into the adapter instead of the insert.
rule cutadapt_trim:
    input:
        r1 = lambda wildcards: PAIRED_SAMPLES[wildcards.sample]
    output:
        trimmed_r1 = TRIMMED_DIR + "/{sample}.cut.R1.fastq"
    params:
        adapter = config["adapters"]["read1"],
        min_length = config["cutadapt"]["min_length"],
        quality_cutoff = config["cutadapt"]["quality_cutoff"],
        trim_5 = config["cutadapt"]["trim_5"],
        trim_3 = config["cutadapt"]["trim_3"],
        trimmed_dir = TRIMMED_DIR
    threads: config["cutadapt"]["cores"]
    resources:
        mem_mb = 24000
    log:
        LOGS_DIR + "/cutadapt/{sample}.log"
    shell:
        """
        module load cutadapt/2.10

        # Create output directory
        mkdir -p {params.trimmed_dir}

        if [ {params.trim_5} -gt 0 ] || [ {params.trim_3} -gt 0 ]; then
            # Two-pass (e.g. NEXTFLEX v3): adapter first, then remove randomers.
            tmp="{output.trimmed_r1}.adapt.tmp"
            cutadapt --quality-cutoff {params.quality_cutoff} \
                --adapter {params.adapter} \
                --cores {threads} \
                --output "$tmp" \
                {input.r1} > {log} 2>&1
            cutadapt -u {params.trim_5} -u -{params.trim_3} \
                --minimum-length {params.min_length} \
                --cores {threads} \
                --output {output.trimmed_r1} \
                "$tmp" >> {log} 2>&1
            rm -f "$tmp"
        else
            # Single-pass (e.g. v4): no randomer end-trimming.
            cutadapt --quality-cutoff {params.quality_cutoff} \
                --adapter {params.adapter} \
                --minimum-length {params.min_length} \
                --cores {threads} \
                --output {output.trimmed_r1} \
                {input.r1} > {log} 2>&1
        fi

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
        OUTPUT_DIR + "/processing_summary.txt"
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

# Rule 6: MultiQC report
rule multiqc:
    input:
        # Depend on trimmed R1 FASTQ outputs
        expand("{trimmed_dir}/{sample}.cut.R1.fastq",
               trimmed_dir=TRIMMED_DIR,
               sample=PAIRED_SAMPLES.keys()),
        # Depend on FastQC zips (raw and trimmed) for MultiQC aggregation
        expand("{fastqc_dir}/{sample}-R1_fastqc.zip",
               fastqc_dir=FASTQC_DIR,
               sample=PAIRED_SAMPLES.keys()),
        # expand("{fastqc_trim_dir}/{sample}.cut.R1_fastqc.zip",
        #        fastqc_trim_dir=FASTQC_DIR + "/trimmed",
        #        sample=PAIRED_SAMPLES.keys()),
        # Include cutadapt logs for MultiQC cutadapt module
        expand("{logs_dir}/cutadapt/{sample}.log",
               logs_dir=LOGS_DIR,
               sample=PAIRED_SAMPLES.keys())
    output:
        report = OUTPUT_DIR + "/multiqc_report.html"
    params:
        output_dir = OUTPUT_DIR
    threads: 2
    resources:
        mem_mb = 4000,
        cpus = 2,
        partition = "standard",
        account = "sbsandme_lab"
    shell:
        """
        rm -f {params.output_dir}/multiqc_report.html || true
        rm -rf {params.output_dir}/multiqc_data || true
        module load singularity/3.11.3
        singularity run /dfs9/ucightf-lab/kstachel/TOOLS/multiqc-1.20.sif \
            multiqc {params.output_dir} -o {params.output_dir}
        module unload singularity/3.11.3
        """

# ---------------------------------------------------------------------------
# Post-CLC reporting
#
# The workflow is interrupted after `multiqc`: the trimmed R1 FASTQs
# (output/trimmed_fastq) are imported into QIAGEN CLC Genomics Workbench (a GUI
# program) and quantified against miRBase. CLC exports per-sample
# "grouped on mature" CSVs into CLC_DIR. The rules below resume from those
# exports and are NOT part of `rule all`; run them once the CLC step is done:
#
#     snakemake report --cores 1
# ---------------------------------------------------------------------------
CLC_DIR = OUTPUT_DIR + "/clc_genomics/v25"
RESULTS_DIR = "results"
CLC_MATURE_CSVS = glob.glob(f"{CLC_DIR}/*grouped on mature*.csv")

# Combine the per-sample CLC miRNA quantifications into long and wide count
# tables. src/parse_miRNA_out.r reads CLC_DIR and writes to RESULTS_DIR (both
# paths are hard-coded in that script).
rule parse_mirna:
    input:
        CLC_MATURE_CSVS
    output:
        long = RESULTS_DIR + "/miRNA_counts_grouped_on_mature_long.csv",
        wide = RESULTS_DIR + "/miRNA_counts_grouped_on_mature.csv"
    log:
        LOGS_DIR + "/parse_mirna.log"
    shell:
        # R + tidyverse/fs/janitor/datapasta are provided by the pixi env
        # (run via `pixi run parse`), so no `module load R` here.
        """
        Rscript src/parse_miRNA_out.r > {log} 2>&1
        """

# Generate the project summary PDF: pipeline description, per-sample metadata
# extracted from the raw FASTQs, and references. Depends on the MultiQC report
# and parsed counts so it only builds once the full project is complete.
rule report:
    input:
        multiqc = OUTPUT_DIR + "/multiqc_report.html",
        counts = RESULTS_DIR + "/miRNA_counts_grouped_on_mature.csv"
    output:
        pdf = OUTPUT_DIR + "/SmallRNA_Project_Report.pdf"
    params:
        fastq_dir = FASTQ_DIR,
        clc_version = "v25",
        library_kit = config.get("library_kit", "Revvity NEXTflex Small RNA v4")
    log:
        LOGS_DIR + "/report.log"
    shell:
        """
        python3 src/generate_report.py \
            --fastq-dir {params.fastq_dir} \
            --output {output.pdf} \
            --clc-version {params.clc_version} \
            --library-kit "{params.library_kit}" > {log} 2>&1
        """
