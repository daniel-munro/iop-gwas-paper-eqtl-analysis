suppressPackageStartupMessages(library(tidyverse))

load_tensorqtl <- function(tensorqtl_out) {
    read_tsv(tensorqtl_out, col_types = "c-----c---------d--") %>%
        # mutate(qval = qvalue::qvalue(pval_beta)$qvalues) %>%
        rename(gene_id = phenotype_id,
               pval = pval_beta) %>%
        separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE,
                 remove = FALSE) %>%
        mutate(chrom = str_replace(chrom, "chr", "") %>% as.integer())
}

load_qtl2 <- function(qtl2_out) {
    read_tsv(qtl2_out, col_types = "cc-dd---------") %>%
        rename(gene_id = gene, variant_id = snp) %>%
        separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE,
                 remove = FALSE) %>%
        mutate(chrom = str_replace(chrom, "chr", "") %>% as.integer())
}

load_afc <- function(afc_out) {
    read_tsv(afc_out, col_types = "cc--d--") %>%
        rename(gene_id = pid,
               variant_id = sid,
               log2_aFC_eQTL = log2_aFC) %>%
        mutate(log2_aFC_eQTL = if_else(is.nan(log2_aFC_eQTL),
                                       as.double(NA),
                                       log2_aFC_eQTL))
}

load_ase <- function(ase_out, min_het_inds = 10) {
    read_tsv(ase_out, col_types = "cc--i---d-------------------") %>%
        filter(var_het_n >= min_het_inds) %>%
        select(gene_id = gene,
               variant_id = var_id,
               log2_aFC_ASE = var_het_afc)
}

load_eqtls <- function(method) {
    d <- if (method == "qtl2") {
        load_qtl2(glue::glue("data/qtl2/eyes.gene_var_pval.tsv.gz"))
    } else {
        load_tensorqtl(glue::glue("data/tensorqtl/{method}.cis_qtl.txt.gz"))
    }
    d %>%
        left_join(load_afc(glue::glue("data/afc/{method}.aFC.txt")),
                  by = c("gene_id", "variant_id")) %>%
        left_join(load_ase(glue::glue("data/afc/{method}.ASE_aFC.txt")),
                  by = c("gene_id", "variant_id"))
}

# eqtls <- tibble(method = c("basic", "basic3", "basic4", "main", "main3", "main4", "main5")) %>%
eqtls <- tibble(method = c("basic4", "main4", "main5")) %>%
    group_by(method) %>%
    summarise({
        # d <- if (method == "qtl2") {
        #     load_qtl2(glue::glue("data/qtl2/eyes.gene_var_pval.tsv.gz"))
        # } else {
        #     load_tensorqtl(glue::glue("data/tensorqtl/{method}.cis_qtl.txt.gz"))
        # }
        # d %>%
        #     left_join(load_afc(glue::glue("data/afc/{method}.aFC.txt")),
        #               by = c("gene_id", "variant_id")) %>%
        #     left_join(load_ase(glue::glue("data/afc/{method}.ASE_aFC.txt")),
        #               by = c("gene_id", "variant_id"))
        load_eqtls(method)
    }, .groups = "drop")
    # filter(qval < 0.05)

write_tsv(eqtls, "data/eqtls.txt")

eqtls %>%
    mutate(across(c("qval", "log2_aFC_eQTL", "log2_aFC_ASE"), signif, digits = 6)) %>%
    filter(method == "basic4") %>%
    select(-method) %>%
    arrange(pval) %>%
    write_tsv("data/eyes.eqtls.txt")

