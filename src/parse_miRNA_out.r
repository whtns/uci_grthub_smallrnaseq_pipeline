#!/usr/bin/env Rscript
library(tidyverse)
library(fs)
library(janitor)
library(datapasta)

long_counts_table <- dir_ls("output/clc_genomics_workbench", glob = "*mature*")  |> 
    map_df(~{
        read_csv(.x) |>
        mutate(sample = path_file(.x) |> str_remove(".cut.R1.*"))
    }) |>
    janitor::clean_names()  |> 
    dplyr::select(name, sample, expression_values)

## write a wide table containing values normalized by hsa-miR-103a-3p
counts_wide <- 
    long_counts_table |> 
    dplyr::select(name, sample, expression_values) |> 
    tidyr::pivot_wider(names_from = "sample", values_from = "expression_values")

write_csv(counts_wide, "results/miRNA_counts_grouped_on_mature0.csv")


## compute per-sample total counts and the reference miRNA (hsa-miR-103a-3p)
mir103a_ref <-
    long_counts_table |> 
    dplyr::filter(name == "hsa-miR-103a-3p") |> 
    dplyr::select(sample, mir103a_counts = expression_values)

## join reference counts back and compute both CPM and CPM normalized to miR-103a
normalized_counts_table <- 
    long_counts_table  |> 
    dplyr::left_join(mir103a_ref, by = "sample") |> 
    group_by(sample) |>
    mutate(total_counts = sum(expression_values)) |>
    ungroup() |>
    mutate(
        normalized_expression = (expression_values / total_counts) * 1e4,
        normalized_by_mir103a = dplyr::case_when(
            is.na(mir103a_counts) ~ NA_real_,
            mir103a_counts == 0 ~ NA_real_,
            TRUE ~ (expression_values / mir103a_counts) * 1e4
        )
    )

## warn if any samples lack the reference miRNA or have zero counts for it
missing_ref <- mir103a_ref |> dplyr::filter(is.na(mir103a_counts))
zero_ref <- mir103a_ref |> dplyr::filter(mir103a_counts == 0)
if (nrow(missing_ref) > 0) {
    message("Warning: the reference miRNA 'hsa-miR-103a-3p' is missing for the following samples: ",
            paste(missing_ref$sample, collapse = ", "))
}
if (nrow(zero_ref) > 0) {
    message("Warning: the reference miRNA 'hsa-miR-103a-3p' has zero counts for the following samples: ",
            paste(zero_ref$sample, collapse = ", "))
}

## write a wide table containing values normalized by hsa-miR-103a-3p
normalized_by_mir103a_wide <- 
    normalized_counts_table |> 
    dplyr::select(name, sample, normalized_by_mir103a) |> 
    tidyr::pivot_wider(names_from = "sample", values_from = "normalized_by_mir103a")

write_csv(normalized_by_mir103a_wide, "results/miRNA_counts_table_normalized_by_hsa-miR-103a-3p.csv")


top_mirnas <- c("hsa-miR-517c-3p", "hsa-miR-517-5p", "hsa-miR-525-3p", "hsa-miR-526b-3p")


ggplot(normalized_by_mir103a_wide  |> 
    dplyr::filter(name %in% top_mirnas) |>
    tidyr::pivot_longer(-name, names_to = "sample", values_to = "counts") ,
    aes(x = name, y = counts + 1, fill = sample)
) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_y_log10() +
    labs(y = "Counts (log10 scale)", x = "miRNA", fill = "sample") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    labs(title = glue::glue("Top miRNA counts normalized by hsa-miR-103a-3p * {1e4}")) +
    NULL


ggsave("results/miRNA_top_counts_normalized_by_hsa-miR-103a-3p.png", width = 8, height = 6)
