signif <- read_tsv("data/Eye.cis_qtl_signif.txt.gz", col_types = "ccidiidddd") |>
    rename(gene_id = phenotype_id) |>
    separate_wider_delim(variant_id, ":", names = c("chrom", "pos"),
             cols_remove = FALSE) |>
    mutate(chrom = str_replace(chrom, "chr", ""),
           pos = as.integer(pos)) |>
    select(-af, -ma_samples, -ma_count, -slope, -slope_se)

tyr <- signif |>
    filter(gene_id == "ENSRNOG00000016421")

angpt2 <- signif |>
    filter(gene_id == "ENSRNOG00000016696")

Rbm12b <- signif |>
    filter(gene_id == "ENSRNOG00000016330")
write_tsv(Rbm12b, "analysis/signif_Rbm12b.tsv")

LOC679087 <- signif |>
    filter(gene_id == "ENSRNOG00000055292")
write_tsv(LOC679087, "analysis/signif_LOC679087.tsv")

Ctsc <- signif |>
    filter(gene_id == "ENSRNOG00000016496")
write_tsv(Ctsc, "analysis/signif_Ctsc.tsv")
