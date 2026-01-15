#!/usr/bin/env Rscript
library(tidyverse)
library(fs)
library(janitor)
library(datapasta)

long_counts_table <- dir_ls("output/clc_genomics", glob = "*mature*")  |> 
    map_df(~{
        read_csv(.x) |>
        mutate(sample = path_file(.x) |> str_remove(".cut.R1.*"))
    }) |>
    janitor::clean_names()  |> 
    dplyr::select(feature_id, name, sample, expression_values)

write_csv(long_counts_table, "results/miRNA_counts_grouped_on_mature_long.csv")

## write a wide table containing values normalized by hsa-miR-103a-3p
counts_wide <- 
    long_counts_table |> 
    dplyr::select(feature_id, name, sample, expression_values) |> 
    tidyr::pivot_wider(names_from = "sample", values_from = "expression_values")

counts_wide[is.na(counts_wide)] <- 0

write_csv(counts_wide, "results/miRNA_counts_grouped_on_mature_wide.csv")
