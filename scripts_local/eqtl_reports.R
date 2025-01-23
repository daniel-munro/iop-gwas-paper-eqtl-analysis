library(tidyverse)

genes <- read_tsv("data/genes.txt", col_types = "cc---ci-----")

alleles <- read_tsv("data/old_unfilt_geno/alleles.txt.gz", col_types = "ccc",
                    col_names = c("variant_id", "ref", "alt"))

afc <- read_tsv("data/Eye.aFC.txt", col_types = "cc--d--") |>
    rename(gene_id = pid,
           variant_id = sid) |>
    mutate(log2_aFC = signif(log2_aFC, 6))

top_assoc <- read_tsv("data/Eye.cis_qtl.txt.gz", col_types = "ci----c---dddd-ddd") |>
    rename(gene_id = phenotype_id) |>
    separate_wider_delim(variant_id, ":", names = c("chrom", "pos"),
                         cols_remove = FALSE) |>
    mutate(chrom = str_replace(chrom, "chr", ""),
           pos = as.integer(pos)) |>
    left_join(afc, by = c("gene_id", "variant_id"), relationship = "one-to-one") |>
    left_join(genes, by = "gene_id", relationship = "one-to-one") |>
    mutate(tss_distance = if_else(strand == "+", pos - tss, tss - pos)) |>
    select(-strand, -tss) |>
    left_join(alleles, by = "variant_id", relationship = "many-to-one") |>
    relocate(gene_name, .after = gene_id) |>
    relocate(tss_distance, .after = af) |>
    relocate(ref, alt, .after = pos)

write_tsv(top_assoc, "data/Eye.top_assoc.txt")

# For paper, top eQTLs
top_assoc |>
    select(-qval) |>
    filter(pval_nominal < pval_nominal_threshold) |>
    write_tsv("data/Eye.top_eQTLs.txt")

eqtls_ind <- read_tsv("data/Eye.cis_independent_qtl.txt.gz",
                      col_types = "ci----c---dddd-cd") |>
    rename(gene_id = phenotype_id) |>
    separate_wider_delim(variant_id, ":", names = c("chrom", "pos"),
             cols_remove = FALSE) |>
    mutate(chrom = str_replace(chrom, "chr", ""),
           pos = as.integer(pos)) |>
    left_join(afc, by = c("gene_id", "variant_id"), relationship = "one-to-one") |>
    left_join(genes, by = "gene_id", relationship = "many-to-one") |>
    mutate(tss_distance = if_else(strand == "+", pos - tss, tss - pos)) |>
    select(-strand, -tss) |>
    left_join(alleles, by = "variant_id", relationship = "many-to-one") |>
    relocate(gene_name, .after = gene_id) |>
    relocate(tss_distance, .after = af) |>
    relocate(ref, alt, .after = pos)

write_tsv(eqtls_ind, "data/Eye.eqtls_indep.txt")
