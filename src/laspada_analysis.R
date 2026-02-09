#!/usr/bin/env Rscript 
# DESeq2 analysis script for miRNA differential expression
#
# USAGE:
#   Rscript laspada_analysis.R [counts_file] [metadata_file] [output_dir]
#
# ARGUMENTS (all optional):
#   counts_file    - CSV file with miRNA counts (rows=miRNAs, cols=samples)
#                    Default: results/miRNA_counts_grouped_on_mature.csv
#   metadata_file  - CSV file with sample metadata (must have sample_name column)
#                    Default: metadata/metadata.csv
#   output_dir     - Directory for output files
#                    Default: output/deseq2
#
# COMPARISONS:
#   Comparisons are defined in metadata/comparisons.yaml
#   Format: list of [group_a, group_b] pairs, e.g.:
#     - [control, treated]
#     - [control, diseased]
#   If file doesn't exist, defaults to: - [control, treated]
#
# OUTPUTS:
#   - dds.rds: DESeq2 dataset object
#   - all_sample_counts.csv: normalized counts for all samples
#   - results_<comparison>.csv: DE results for each comparison
#   - results_<comparison>.rds: DESeqResults object for each comparison
#   - results_all_contrasts.csv: combined results from all comparisons
#   - results/pca_plots.pdf: PCA plots
#
# EXAMPLES:
#   # Run with all defaults
#   Rscript laspada_analysis.R
#
#   # Custom counts file only
#   Rscript laspada_analysis.R my_counts.csv
#
#   # Custom counts and metadata
#   Rscript laspada_analysis.R my_counts.csv my_metadata.csv
#
#   # All custom
#   Rscript laspada_analysis.R my_counts.csv my_metadata.csv custom_output/

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(glue)
  library(AnnotationDbi)
  library(yaml)
  library(iSEEde)
  library(purrr)
})

args <- commandArgs(trailingOnly=TRUE)
default_counts <- "results/miRNA_counts_grouped_on_mature.csv"
## Defaults (allow running without providing all args)
default_meta <- "metadata/metadata.csv"
default_out <- "output/deseq2"

# Use provided counts file if given, otherwise fall back to default
counts_file <- default_counts
if(length(args) >= 1 && nzchar(args[1])) counts_file <- args[1]

# meta and out dir: use provided values if present, otherwise use defaults
meta_file <- default_meta
if(length(args) >= 2 && nzchar(args[2])) meta_file <- args[2]

out_dir <- default_out
if(length(args) >= 3 && nzchar(args[3])) out_dir <- args[3]

# Optional arguments for condition column name and specific groups
condition_col <- args[4]  # Optional: condition column name in metadata
if(is.na(condition_col) || !nzchar(condition_col)) condition_col <- "condition"

dir.create(out_dir, showWarnings=FALSE)

counts <- read.csv(counts_file, header=TRUE, check.names=FALSE)

# Use name column as row identifiers (these are the miRNA names)
if("name" %in% colnames(counts)){
  count_matrix <- counts[, setdiff(colnames(counts), c("feature_id", "name")), drop=FALSE]
  row_ids <- counts$name
  # Make row names unique by appending sequence number if duplicates exist
  row_ids <- make.unique(row_ids, sep=".")
  rownames(count_matrix) <- row_ids
} else if("feature_id" %in% colnames(counts)){
  count_matrix <- counts[, -which(colnames(counts) %in% c("feature_id")), drop=FALSE]
  row_ids <- counts$feature_id
  row_ids <- make.unique(row_ids, sep=".")
  rownames(count_matrix) <- row_ids
} else {
  # Use first column as row names
  count_matrix <- counts[, -1, drop=FALSE]
  row_ids <- make.unique(as.character(counts[[1]]), sep=".")
  rownames(count_matrix) <- row_ids
}

# Ensure numeric matrix
count_matrix <- as.matrix(count_matrix)
mode(count_matrix) <- "integer" 

meta <- read.csv(meta_file)  |> 
janitor::clean_names() |>
dplyr::mutate(sample_name = str_replace_all(sample_name, "-", ".")) |>
tibble::column_to_rownames("sample_name") |>
  identity()

# Ensure the condition column exists
if(!condition_col %in% colnames(meta)){
  stop(glue::glue("Condition column '{condition_col}' not found in metadata. Available columns: {paste(colnames(meta), collapse=', ')}"))
}

# Convert condition column to factor if it isn't already
meta[[condition_col]] <- factor(meta[[condition_col]])

colnames(count_matrix) <- rownames(meta)

# Write normalized counts before converting to matrix
count_df <- as.data.frame(count_matrix) |> 
  tibble::rownames_to_column("miRNA") |>
  readr::write_csv(file.path(out_dir, "all_sample_counts.csv"))

print(rownames(meta))
print(colnames(count_matrix))
# Ensure sample names match
count_matrix <- count_matrix[, rownames(meta)]

# Create design formula dynamically
design_formula <- as.formula(paste0("~", condition_col))
dds <- DESeqDataSetFromMatrix(countData=count_matrix, 
  colData=meta, 
  design=design_formula)

dds <- DESeq(dds)


# 6. Quality Control and Visualization
# ----------------------------------------------------------- #
# Optional: Transform data for visualization
# Try VST first, fall back to rlog if there are too few genes
transformed_data <- tryCatch({
  vst(dds, blind = FALSE) # Variance stabilizing transformation
}, error = function(e){
  message(glue::glue("VST failed: {e$message}. Using rlog instead."))
  rlog(dds, blind = FALSE)
})

plot_var <- c(condition_col) |> 
set_names()

pca_plots <- map(plot_var, ~{plotPCA(transformed_data, intgroup = .x) + labs(title = glue("PCA - {.x}"))})

pdf(file.path(out_dir, "pca_plots.pdf"))
print(pca_plots)
dev.off()

comparisons_yaml <- file.path(dirname(meta_file), "comparisons.yaml")

# If a comparisons YAML exists, read it; otherwise create a default and write it.
if(file.exists(comparisons_yaml)){
  comparisons <- tryCatch({
    yaml::read_yaml(comparisons_yaml)
  }, error = function(e){
    warning(glue::glue("Failed to read {comparisons_yaml}: {e$message}. Falling back to defaults."))
    NULL
  })
  # validate structure
  if(is.null(comparisons) || !is.list(comparisons) || length(comparisons) == 0){
    warning(glue::glue("Invalid comparisons in {comparisons_yaml}; using defaults."))
    comparisons <- list(c("control", "treated"))
    try({ yaml::write_yaml(comparisons, comparisons_yaml) }, silent = TRUE)
  }
} else {
  comparisons <- list(c("control", "treated"))
  dir.create(dirname(meta_file), showWarnings = FALSE, recursive = TRUE)
  tryCatch({
    yaml::write_yaml(comparisons, comparisons_yaml)
    message(glue::glue("Wrote comparisons to {comparisons_yaml}"))
  }, error = function(e){
    warning(glue::glue("Failed to write comparisons YAML: {e$message}"))
  })
}

# Iterate over defined comparisons, extract results for each, and save RDS + CSV
safe_filename <- function(x){
  x <- stringr::str_replace_all(x, "\\+", "plus")
  x <- stringr::str_replace_all(x, "[^A-Za-z0-9_-]", "_")
  tolower(x)
}

all_res_list <- list()
for(cmp in comparisons){
  group_a <- cmp[1]
  group_b <- cmp[2]

  message(glue::glue("Extracting results: {group_a} vs {group_b}"))

  # results() uses contrast = c("factorName","level1","level2")
  res_i <- tryCatch(
    results(dds, contrast = c(condition_col, group_a, group_b)),
    error = function(e){
      warning(glue::glue("Failed to get results for {group_a} vs {group_b}: {e$message}"))
      NULL
    }
  )

  if(is.null(res_i)) next
  # filter out rows where pvalue is NA
  res_df0 <- as.data.frame(res_i)
  if("pvalue" %in% colnames(res_df0)){
    res_df <- res_df0[!is.na(res_df0$pvalue), , drop = FALSE]
  } else {
    res_df <- res_df0
  }

  if(nrow(res_df) == 0){
    message(glue::glue("No rows with non-NA pvalue for {group_a} vs {group_b}; skipping."))
    next
  }

  nm <- glue::glue("{group_a}_vs_{group_b}") %>% safe_filename()
  rds_path <- file.path(out_dir, glue::glue("results_{nm}.rds"))
  csv_path <- file.path(out_dir, glue::glue("results_{nm}.csv"))

  # save filtered results: save RDS as DESeqResults subset and CSV as data.frame
  keep_idx <- rownames(res_df)
  res_i_filt <- tryCatch(res_i[keep_idx, ], error = function(e) res_i)
  saveRDS(res_i_filt, file = rds_path)

  # store for combined table (add contrast column)
  df_i <- res_df
  df_i$miRNA <- rownames(df_i)
  df_i$contrast <- glue::glue("{group_a}_vs_{group_b}")

  # reorder columns: miRNA first, then everything else
  other_cols <- setdiff(colnames(df_i), c("miRNA"))
  df_i <- df_i[, c("miRNA", other_cols), drop = FALSE]

  # sort by increasing padj (NA last) if padj column exists
  if("padj" %in% colnames(df_i)){
    df_i <- df_i[order(df_i$padj, na.last = TRUE), , drop = FALSE]
  }

  write.csv(df_i, file = csv_path, row.names = FALSE)
  all_res_list[[nm]] <- df_i
}

# Save combined results if any
if(length(all_res_list) > 0){
  combined <- dplyr::bind_rows(all_res_list)
  if("padj" %in% colnames(combined)){
    combined <- combined[order(combined$padj, na.last = TRUE), , drop = FALSE]
  }
  write.csv(combined, file = file.path(out_dir, "results_all_contrasts.csv"), row.names = FALSE)
}

# Save the full DESeq dataset for downstream use
saveRDS(dds, file=file.path(out_dir, "dds.rds"))
