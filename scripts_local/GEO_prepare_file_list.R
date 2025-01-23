library(tidyverse)

fastq <- read_tsv("data/fastq_map.txt", col_types = "ccc",
                  col_names = c("fastq1", "fastq2", "sample"))

fastq_removed <- read_tsv("data/fastq_map.pre_qc.txt", col_types = "ccc",
                          col_names = c("fastq1", "fastq2", "sample")) |>
    filter(!(fastq1 %in% fastq$fastq1)) |>
    mutate(sample = str_c("Removed_", as.integer(fct_inorder(sample))))

## File list and description per sample to paste into metadata spreadsheet

df <- fastq |>
    bind_rows(fastq_removed) |>
    pivot_longer(c(fastq1, fastq2), names_to = "end", values_to = "fastq") |>
    mutate(fastq = str_replace(fastq, "/s", "_s")) |>
    mutate(file_n = str_c("file", 1:n()), .by = sample) |>
    pivot_wider(id_cols = sample, names_from = file_n, values_from = fastq) |>
    mutate(
        original_id = str_split_fixed(file1, "_", 5)[, 4],
        description = str_glue("Original label before genotype mismatch fix: {original_id}"),
        description = if_else(sample == original_id, "", as.character(description)),
        description = if_else(str_sub(sample, 1, 7) == "Removed",
                              "Removed: no genotype match found",
                              description)
    )

write_tsv(df, "qc/GEO/sample_files.txt", col_names = FALSE, na = "")

## List of paired-end file pairs to paste into metadata spreadsheet

paired <- fastq |>
    bind_rows(fastq_removed) |>
    mutate(fastq1 = str_replace(fastq1, "/s", "_s"),
           fastq2 = str_replace(fastq2, "/s", "_s"))

write_tsv(paired, "qc/GEO/paired_files.txt", col_names = FALSE)

## Commands to collapse FASTQ files into one folder (as hardlink copies)

ln_commands <- fastq |>
    bind_rows(fastq_removed) |>
    pivot_longer(c(fastq1, fastq2), names_to = "end", values_to = "fastq") |>
    mutate(fastq_new = str_replace(fastq, "/s", "_s"),
           command = str_glue("ln {fastq} GEO/fastq/{fastq_new}"))

ln_commands |>
    select(command) |>
    write_tsv("qc/GEO/ln_commands.txt", col_names = FALSE)
