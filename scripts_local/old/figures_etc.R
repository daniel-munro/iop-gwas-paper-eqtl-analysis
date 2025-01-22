library(tidyverse)
library(patchwork)

genes <- read_tsv("data/reference/GENES_RAT.txt",
                  col_types = cols(ENSEMBL_ID = "c", SYMBOL = "c", .default = "-"),
                  skip = 83) %>%
    rename(gene_id = ENSEMBL_ID,
           gene = SYMBOL) %>%
    separate_rows(gene_id, sep=";") %>%
    filter(!str_detect(gene, "^LOC\\d+$")) %>%
    group_by(gene_id) %>%
    slice(1) %>%
    ungroup()

signif <- read_tsv("data/tensorqtl/main5.cis_qtl_signif.txt.gz",
                   col_types = "ccidiid--d") %>%
    rename(gene_id = phenotype_id) %>%
    separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE,
             remove = FALSE) %>%
    mutate(chrom = str_replace(chrom, "chr", "") %>% as.integer(),
           pos = pos / 1e6)

perm <- read_tsv("data/tensorqtl/main5.cis_qtl.txt.gz",
                 col_types = "ci----cd----d---ddd") %>%
    rename(gene_id = phenotype_id) %>%
    separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE,
             remove = FALSE) %>%
    mutate(chrom = str_replace(chrom, "chr", "") %>% as.integer(),
           tss = (pos - tss_distance) / 1e6,
           pos = pos / 1e6) %>%
    left_join(genes, by = "gene_id")

iop <- signif %>%
    filter(chrom == 16,
           pos > 76.543372 - 1,
           pos < 76.543372 + 1)

chr16_all <- read_tsv("data/tensorqtl/main5.cis_qtl_pairs.16.txt.gz",
                      col_types = "ccidiiddd") %>%
    rename(gene_id = phenotype_id) %>%
    separate(variant_id, c("chrom", "pos"), sep=":", convert = TRUE,
             remove = FALSE) %>%
    mutate(tss = (pos - tss_distance) / 1e6,
           pos = pos / 1e6) %>%
    left_join(genes, by = "gene_id")

iop_all <- chr16_all %>%
    filter(pos > 76.543372 - 5.1,
           pos < 76.543372 + 5.1)

# tss <- iop_all %>%
#     group_by(gene_id, gene) %>%
#     slice(1) %>%
#     summarise(tss = pos - tss_distance, .groups = "drop") %>%
#     left_join(select(perm, gene_id, pval_nominal_threshold), by = "gene_id") %>%
#     mutate(mlog10_threshold = -log10(pval_nominal_threshold))

iop_all %>%
    # slice_sample(n = 1000) %>%
    mutate(`-log10(pval_nominal)` = -log10(pval_nominal)) %>%
    ggplot(aes(x = pos, y = `-log10(pval_nominal)`, color = gene)) +
    geom_vline(xintercept = 76.543372, color = "#555555") +
    geom_vline(aes(xintercept = tss, color = gene), data = tss) +
    # geom_hline(aes(yintercept = mlog10_threshold, color = gene_id),
    #            data = filter(perm, gene_id %in% tss$gene_id),
    #            alpha = 0.5) +
    geom_point(aes(x = tss, y = mlog10_threshold, color = gene),
               data = tss,
               shape = 4) +
    geom_point(size = 0.5) +
    coord_cartesian(xlim = c(76.543372 - 1, 76.543372 + 1))

ggsave("figures/IOP_QTL_region.png", width = 8, height = 4)

######################################
## Top assoc for genes near IOP QTL ##
######################################

## Use function below instead
# p1 <- perm %>%
#     mutate(pos = tss) %>%
#     filter(chrom == 16,
#            pos > 76.543372 - 3.1,
#            pos < 76.543372 + 3.1) %>%
#     mutate(log10_pval = -log10(pval_nominal),
#            log10_threshold = -log10(pval_nominal_threshold),
#            signif = pval_nominal < pval_nominal_threshold,
#            gene = if_else(pval_nominal < 0.001 | gene == "Angpt2", gene, "")) %>%
#     ggplot(aes(x = pos, y = log10_pval, label = gene)) +
#     geom_segment(aes(x = pos, y = log10_pval,
#                      xend = pos, yend = log10_threshold),
#                  size = 0.3, alpha = 0.5, color = "blue") +
#     geom_vline(xintercept = 76.543372, color = "red", lty = 3) +
#     geom_point(size = 1) +
#     ggrepel::geom_text_repel(force = 10, max.iter = 1e4) +
#     annotate("text", x = 76.6, y = 0.5, color = "red", label = "IOP QTL",
#              hjust = 0) +
#     expand_limits(y = 0) +
#     coord_cartesian(xlim = c(73.54, 79.54), expand = FALSE) +
#     xlab("chr16 position (Mbp)") +
#     ylab(expression(-log[10]("p-value"))) +
#     theme_minimal()
# p1
# 
# ggsave("figures/IOP_QTL_genes.png", width = 6, height = 5)

#########################
## Angpt2 associations ##
#########################

angpt2 <- read_tsv("inspect_qtl/ENSRNOG00000016696.4mbp.cis_qtl_pairs.txt",
                   col_types = "ccidccddd") %>%
    rename(gene_id = phenotype_id) %>%
    separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE,
             remove = FALSE) %>%
    mutate(chrom = str_replace(chrom, "chr", "") %>% as.integer(),
           pos = pos / 1e6)

angpt2_tss <- angpt2$pos[1] - angpt2$tss_distance[1] / 1e6
angpt2_thr <- perm %>%
    filter(gene == "Angpt2") %>%
    pull(pval_nominal_threshold)

p2 <- angpt2 %>%
    mutate(log10_pval = -log10(pval_nominal)) %>%
    ggplot(aes(x = pos, y = log10_pval)) +
    geom_vline(xintercept = angpt2_tss, color = "#777777") +
    geom_hline(yintercept = -log10(angpt2_thr), color = "blue", lty = 2) +
    geom_point(size = 0.3) +
    geom_vline(xintercept = 76.543372, color = "red", lty = 3) +
    xlab("chr16 position (Mbp)") +
    ylab(expression(-log[10]("p-value"))) +
    theme_minimal() +
    annotate("text", x = 75.92, y = 3, color = "#666666", label = "Angpt2\nTSS",
             hjust = 1) +
    annotate("text", x = 76.6, y = 3, color = "red", label = "IOP QTL",
             hjust = 0) +
    annotate("text", x = 73.9, y = 3.95, color = "blue", label = "Angpt2 eQTL threshold",
             hjust = 0) +
    coord_cartesian(xlim = c(73.54, 79.54), expand = FALSE)
p2

ggsave("figures/Angpt2.png")

p1 / p2 + plot_layout(heights = c(2, 1))
ggsave("figures/eQTL_figure.png", width = 7, height = 6)


#####################################
## Like panel A but for 3 new QTLs ##
#####################################

plot_qtl_region <- function(chr, posn) {
    perm %>%
        mutate(pos = tss) %>%
        filter(chrom == chr,
               pos > posn - 3.1,
               pos < posn + 3.1) %>%
        mutate(log10_pval = -log10(pval_nominal),
               log10_threshold = -log10(pval_nominal_threshold),
               signif = pval_nominal < pval_nominal_threshold,
               gene = if_else(pval_nominal < 0.001 | gene == "Angpt2", gene, "")) %>%
        ggplot(aes(x = pos, y = log10_pval, label = gene)) +
        geom_segment(aes(x = pos, y = log10_pval,
                         xend = pos, yend = log10_threshold),
                     size = 0.3, alpha = 0.5, color = "blue") +
        geom_vline(xintercept = posn, color = "red", lty = 3) +
        geom_point(size = 1) +
        ggrepel::geom_text_repel(force = 10, max.iter = 1e4) +
        annotate("text", x = posn + 0.1, y = 0.5, color = "red", label = "IOP QTL",
                 hjust = 0) +
        expand_limits(y = 0) +
        coord_cartesian(xlim = c(posn - 3, posn + 3), expand = FALSE) +
        xlab(str_glue("chr{chr} position (Mbp)")) +
        ylab(expression(-log[10]("p-value"))) +
        theme_minimal()
}
plot_qtl_region(16, 76.543372)
plot_qtl_region(16, 76.528844)
p1 <- plot_qtl_region(1, 151.709252)
p2 <- plot_qtl_region(5, 23.768259)
p1 / p2

ggsave("figures/new_IOP_QTLs.png", width = 6, height = 5)
