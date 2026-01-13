#!/usr/bin/env Rscript
library(tidyverse)
library(fs)
library(janitor)
library(datapasta)

p066 <- read_csv("output/clc_genomics_workbench/xR074-L7-G6-P066-ACCTATAT-TCAGCTGG.csv") |> 
dplyr::mutate(sample = "xR074-L7-G6-P066-ACCTATAT-TCAGCTGG")  |> 
    janitor::clean_names() |> 
    dplyr::select(name, sample, expression_values)  |> 
    identity()

long_counts_table <- read_csv("output/clc_genomics_workbench/grouped_on_mature.csv") |> 
    janitor::clean_names()  |> 
    dplyr::select(name, sample, expression_values)  |> 
    dplyr::bind_rows(p066) |>
    dplyr::mutate(sample = stringr::str_remove(sample, ".cut.*"))  |> 
    dplyr::arrange(sample)

## write a wide table containing values normalized by hsa-miR-103a-3p
counts_wide <- 
    long_counts_table |> 
    dplyr::select(name, sample, expression_values) |> 
    tidyr::pivot_wider(names_from = "sample", values_from = "expression_values")

write_csv(counts_wide, "results/miRNA_counts_grouped_on_mature.csv")

