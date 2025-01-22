signif <- read_tsv("data/Eye.cis_qtl_signif.txt.gz", col_types = "ccidiidddd") |>
    rename(gene_id = phenotype_id) |>
    separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE,
             remove = FALSE) |>
    mutate(chrom = str_replace(chrom, "chr", "")) |>
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

# Use below if you want more than cis-window SNPs, not really valid though

# # Extract eSNPs for Tyr and Angpt2 for Apurva to check LD with IOP QTL or albinism causal SNP.
# # Note that some are not really cis-eSNPs as they are outside the tested cis-window.
# 
# library(tidyverse)
# 
# thr <- read_tsv("data/Eye.cis_qtl.txt.gz", col_types = "c----------------d")
# 
# tyr <- read_tsv("data/inspect_qtl/ENSRNOG00000016421.4mbp.cis_qtl_pairs.txt.gz",
#                 col_types = "cci---d--") |>
#     left_join(thr, by = "phenotype_id") |>
#     rename(gene_id = phenotype_id) |>
#     filter(pval_nominal < pval_nominal_threshold) |>
#     separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE, remove = FALSE) |>
#     mutate(chrom = str_replace(chrom, "chr", ""),
#            pval_nominal = signif(pval_nominal, 6))
# 
# angpt2 <- read_tsv("data/inspect_qtl/ENSRNOG00000016696.4mbp.cis_qtl_pairs.txt.gz",
#                    col_types = "cci---d--") |>
#     left_join(thr, by = "phenotype_id") |>
#     rename(gene_id = phenotype_id) |>
#     filter(pval_nominal < pval_nominal_threshold) |>
#     separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE, remove = FALSE) |>
#     mutate(chrom = str_replace(chrom, "chr", ""),
#            pval_nominal = signif(pval_nominal, 6))
# 
# write_tsv(tyr, "analysis/signif_Tyr.tsv")
# write_tsv(angpt2, "analysis/signif_Angpt2.tsv")
