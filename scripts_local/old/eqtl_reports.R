library(tidyverse)

afc <- read_tsv("data/afc/main5.aFC.txt", col_types = "cc--d--") %>%
    rename(gene_id = pid,
           variant_id = sid) %>%
    mutate(log2_aFC = if_else(is.nan(log2_aFC), as.double(NA), signif(log2_aFC, 6)))

perm <- read_tsv("data/tensorqtl/main5.cis_qtl.txt.gz",
                 col_types = "ci----c-----d---ddd") %>%
    rename(gene_id = phenotype_id) %>%
    separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE,
             remove = FALSE) %>%
    mutate(chrom = str_replace(chrom, "chr", "") %>% as.integer()) %>%
    arrange(pval_beta)

signif <- read_tsv("data/tensorqtl/main5.cis_qtl_signif.txt.gz",
                   col_types = "ccidiid--d") %>%
    rename(gene_id = phenotype_id) %>%
    separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE,
             remove = FALSE) %>%
    mutate(chrom = str_replace(chrom, "chr", "") %>% as.integer()) %>%
    left_join(afc, by = c("gene_id", "variant_id"))

write_tsv(perm, "data/eyes.top_assoc_per_gene.txt")

write_tsv(signif, "data/eyes.all_signif_eqtls.txt")
