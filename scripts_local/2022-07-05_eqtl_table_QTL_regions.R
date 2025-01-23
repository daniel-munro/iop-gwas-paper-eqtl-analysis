# From Oksana: Could you please give me the values for each gene in 3 QTL regions for cis QTLs?
# gene name
# top SNP fro eQTL
# p-value
# significance threshold
# Is there anything else to include in the table?

library(tidyverse)

iop <- c("1" = 151709252, "5" = 23768259, "16" = 76528844)
iop_min <- c("1" = 150500000, "5" = 23000000, "16" = 75500000)  # Figure limits
iop_max <- c("1" = 153500000, "5" = 26000000, "16" = 78000000)  # Figure limits

## Use RGD genes, same info used for figures
genes_all <- read_tsv("data/reference/GENES_RAT.txt",
                  col_types = cols(ENSEMBL_ID = "c", SYMBOL = "c", CHROMOSOME_6.0 = "c",
                                   START_POS_6.0 = "c", STOP_POS_6.0 = "c",
                                   STRAND_6.0 = "c", .default = "-"),
                  skip = 83) |>
    rename(gene_id = ENSEMBL_ID,
           gene_name = SYMBOL,
           chrom = CHROMOSOME_6.0,
           start = START_POS_6.0,
           stop = STOP_POS_6.0,
           strand = STRAND_6.0) |>
    separate_rows(gene_id, sep = ";") |>
    filter(!is.na(gene_id),
           !is.na(start)) |>
    group_by(gene_id) |>
    slice(1) |>
    ungroup() |>
    # For genes with multiple locations, keep first:
    mutate(chrom = str_split(chrom, ";", simplify = TRUE)[, 1],
           strand = str_split(strand, ";", simplify = TRUE)[, 1],
           start = str_split(start, ";", simplify = TRUE)[, 1] |> as.integer(),
           stop = str_split(stop, ";", simplify = TRUE)[, 1] |> as.integer(),
           tss = if_else(strand == "+", start, stop))

genes <- genes_all |>
    filter(
        (chrom == "1" & tss > iop_min["1"] & tss < iop_max["1"]) |
            (chrom == "5" & tss > iop_min["5"] & tss < iop_max["5"]) |
            (chrom == "16" & tss > iop_min["16"] & tss < iop_max["16"])
    )
gene_names <- genes |>
    select(gene_id, gene_name) |>
    deframe()
gene_tss <- genes |>
    select(gene_id, tss) |>
    deframe()

eqtls <- read_tsv("data/Eye.top_assoc.txt", col_types = "ccicciccdiddddddd") |>
    select(gene_id, chrom, variant_id, ref, alt, log2_aFC,
           pval_nominal, pval_nominal_threshold) |>
    filter(gene_id %in% genes$gene_id) |>
    mutate(gene_name = gene_names[gene_id], .after = gene_id) |> # Get RGD gene names to match figures
    mutate(tss = gene_tss[gene_id], .after = chrom) |>
    # In rare case genes with multiple tested gene IDs, keep most significant
    group_by(gene_name) |>
    slice_min(pval_nominal, n = 1) |>
    ungroup() |>
    arrange(as.integer(chrom), tss)

write_tsv(eqtls, "data/analysis/eqtl_table_QTL_regions.txt")
