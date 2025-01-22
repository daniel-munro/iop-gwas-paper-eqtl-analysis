library(tidyverse)
library(ggrepel)

# chr <- 16
# pos <- 76528844
# chr <- 1
# pos <- 151709252
chr <- 5
pos <- 23768259

trans <- read_tsv(str_glue("inspect_qtl/chr{chr}_{pos}.trans_qtl_pairs.txt"),
                  col_types = "ccdddd") %>%
    arrange(pval) %>%
    mutate(qval = qvalue::qvalue(pval)$qvalues)
cis <- read_tsv(str_glue("inspect_qtl/chr{chr}_{pos}.cis_qtl_pairs.txt"),
                col_types = "ccidiiddd") %>%
    arrange(pval_nominal)

write_tsv(trans, str_glue("inspect_qtl/chr{chr}_{pos}.trans_pairs.txt"))
write_tsv(cis, str_glue("inspect_qtl/chr{chr}_{pos}.cis_pairs.txt"))

cis %>%
    mutate(log10_pval = -log10(pval_nominal),
           tss = (pos - tss_distance) / 1e6,
           phenotype_id = if_else(pval_nominal < 0.05, phenotype_id, "")) %>%
    ggplot(aes(x = tss, y = log10_pval, label = phenotype_id)) +
    geom_point() +
    xlab("Gene TSS (Mbp)") +
    ylab("-log10(p-value)") +
    geom_text_repel()

cis %>%
    select(phenotype_id, tss_distance, pval_nominal)

trans %>%
    mutate(qval = qvalue::qvalue(pval)$qvalues) %>%
    arrange(pval) %>%
    select(variant = variant_id,
           gene = phenotype_id,
           pval, qval)

trans %>%
    ggplot(aes(x = pval)) +
    geom_histogram(bins = 30)
