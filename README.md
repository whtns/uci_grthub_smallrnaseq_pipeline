# Small RNA Sequencing Processing Workflow

This Snakemake workflow processes small RNA sequencing data from Biooscientific small RNA v4 kit, performing quality control and adapter trimming as specified in the original shell scripts.

## Overview

The workflow performs the following steps:
1. **FastQC analysis** on raw FASTQ files
2. **Adapter trimming** using Cutadapt with small RNA v4 kit specific adapters
3. **FastQC analysis** on trimmed files
4. **Summary report** generation

![](rulegraph.png}

## Files Structure

```
├── Snakefile                 # Main workflow file
├── config.yaml              # Configuration parameters
├── cluster.yaml             # SLURM cluster configuration
├── FastqFiles/              # Input FASTQ files
├── fastQC/                  # FastQC outputs
│   ├── html/               # Raw FastQC HTML reports
│   └── trimmed_html/       # Trimmed FastQC HTML reports
├── trimmed_fastq/          # Cutadapt trimmed FASTQ files
└── logs/                   # Log files
```

## Requirements

- Snakemake (>=6.0)
- FastQC (0.11.9)
- Cutadapt (2.10)
- SLURM job scheduler (for cluster execution)

## Usage

### Local Execution

```bash
# Dry run to check workflow
snakemake -n

# Run with 4 cores
snakemake --cores 4

# Run only FastQC analysis
snakemake fastqc_only --cores 4

# Run only adapter trimming
snakemake cutadapt_only --cores 4
```

### SLURM Cluster Execution

```bash
# Create logs directory
mkdir -p logs

# Submit jobs to SLURM cluster
snakemake --cluster-config cluster.yaml \
    --cluster "sbatch --account={cluster.account} --partition={cluster.partition} --nodes={cluster.nodes} --ntasks={cluster.ntasks} --cpus-per-task={cluster.cpus} --mem={cluster.mem} --time={cluster.time} --output={cluster.output} --error={cluster.error}" \
    --jobs 10

# Alternative using SLURM profile (if configured)
snakemake --profile slurm --jobs 10
```

### Workflow Options

```bash
# Generate workflow visualization
snakemake --dag | dot -Tpng > workflow_dag.png

# Generate rule graph
snakemake --rulegraph | dot -Tpng > rules.png

# Generate summary report
snakemake create_summary
```

## Input Files

The workflow automatically discovers FASTQ files matching the pattern:
- `FastqFiles/hts.igb.uci.edu/rob25092453/*Sequences.txt.gz`
- Expects paired-end data (READ1 and READ2 files)

Example input files:
- `xR064-L2-G4-P18-AAGGCCTG-TCTCAATT-READ1-Sequences.txt.gz`
- `xR064-L2-G4-P18-AAGGCCTG-TCTCAATT-READ2-Sequences.txt.gz`

## Output Files

### FastQC Reports
- **Raw data**: `fastQC/html/{sample}-READ{1,2}-Sequences_fastqc.html`
- **Trimmed data**: `fastQC/trimmed_html/{sample}.cut.R{1,2}_fastqc.html`

### Trimmed FASTQ Files
- `trimmed_fastq/{sample}.cut.R1.fastq`
- `trimmed_fastq/{sample}.cut.R2.fastq`

### Log Files
- `logs/fastqc_raw/{sample}.log`
- `logs/cutadapt/{sample}.log`
- `logs/fastqc_trimmed/{sample}.log`

## Configuration

Edit `config.yaml` to modify:
- Adapter sequences
- Quality cutoff values
- Minimum read length
- Resource allocation
- Directory paths

## Adapter Information

The workflow uses Biooscientific small RNA v4 kit adapters:
- **Read 1 adapter**: `TGGAATTCTCGGGTGCCAAGG`
- **Read 2 adapter**: `AGATCGGAAGAGCGTCGTGTAGGGAAAGA`

## Quality Control Parameters

- **Quality cutoff**: 20 (removes bases with quality < 20 from 3' end)
- **Minimum length**: 16 bases (excludes reads shorter than 16 bp)
- **Pair adapter mode**: Enabled for proper paired-end trimming

## Troubleshooting

### Common Issues

1. **Module loading errors**: Ensure FastQC and Cutadapt modules are available
2. **File path issues**: Check that input FASTQ files exist in the expected location
3. **Permission errors**: Ensure write permissions for output directories
4. **Memory issues**: Increase memory allocation in `cluster.yaml` if needed

### Checking Results

```bash
# Check number of processed samples
ls trimmed_fastq/*.fastq | wc -l

# View cutadapt statistics
grep -A 20 "Total read pairs processed" logs/cutadapt/*.log

# Check FastQC results
open fastQC/html/*.html
open fastQC/trimmed_html/*.html
```

## Notes

- The workflow is designed for small RNA sequencing data
- Adapter trimming is specifically configured for Biooscientific small RNA v4 kit
- Both reads are processed, but according to the original scripts, analysis typically proceeds with Read 1 only due to small insert sizes
- All intermediate files are preserved for quality control purposes

## References

- [Snakemake Documentation](https://snakemake.readthedocs.io/)
- [FastQC Documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Cutadapt Documentation](https://cutadapt.readthedocs.io/)
- [Biooscientific small RNA v4 Kit](https://www.biosearchtech.com/)# uci_grthub_smallrnaseq_pipeline
