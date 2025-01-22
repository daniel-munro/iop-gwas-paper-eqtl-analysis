library(tidyverse)
library(patchwork)

genes <- read_tsv("data/reference/GENES_RAT.txt",
                  col_types = cols(ENSEMBL_ID = "c", SYMBOL = "c", .default = "-"),
                  skip = 83) |>
    rename(gene_id = ENSEMBL_ID,
           gene = SYMBOL) |>
    separate_rows(gene_id, sep=";") |>
    # filter(!str_detect(gene, "^LOC\\d+$")) |>
    group_by(gene_id) |>
    slice(1) |>
    ungroup()

perm <- read_tsv("data/tensorqtl/main5.cis_qtl.txt.gz",
                 col_types = "ci----cd----d---ddd") |>
    rename(gene_id = phenotype_id) |>
    separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE,
             remove = FALSE) |>
    mutate(chrom = str_replace(chrom, "chr", "") |> as.integer(),
           tss = (pos - tss_distance) / 1e6,
           pos = pos / 1e6) |>
    left_join(genes, by = "gene_id")


#########################
## Angpt2 associations ##
#########################

angpt2 <- read_tsv("inspect_qtl/ENSRNOG00000016696.4mbp.cis_qtl_pairs.txt",
                   col_types = "ccidccddd") |>
    rename(gene_id = phenotype_id) |>
    separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE,
             remove = FALSE) |>
    mutate(chrom = str_replace(chrom, "chr", "") |> as.integer(),
           pos = pos / 1e6)

angpt2_tss <- angpt2$pos[1] - angpt2$tss_distance[1] / 1e6
angpt2_thr <- perm |>
    filter(gene == "Angpt2") |>
    pull(pval_nominal_threshold)

p1 <- angpt2 |>
    filter(!is.na(pval_nominal)) |>
    mutate(log10_pval = -log10(pval_nominal)) |>
    ggplot(aes(x = pos, y = log10_pval)) +
    geom_vline(xintercept = angpt2_tss, color = "#777777") +
    geom_hline(yintercept = -log10(angpt2_thr), color = "blue", lty = 2) +
    geom_point(size = 0.3) +
    geom_vline(xintercept = 76.528844, color = "red", lty = 3) +
    xlab("chr16 position (Mbp)") +
    ylab(expression(-log[10]*"p-value")) +
    theme_minimal() +
    annotate("text", x = 75.93, y = 3, color = "#666666", label = "Angpt2\nTSS",
             hjust = 1) +
    annotate("text", x = 76.578844, y = 3, color = "red", label = "IOP QTL",
             hjust = 0) +
    annotate("text", x = 74.7, y = 3.95, color = "blue", label = "Angpt2 eQTL threshold",
             hjust = 0) +
    coord_cartesian(xlim = c(74.528844, 78.528844), expand = FALSE)
p1

# ggsave("figures/Angpt2.png")

######################
## Tyr associations ##
######################

tyr <- read_tsv("inspect_qtl/ENSRNOG00000016421.4mbp.cis_qtl_pairs.txt",
                col_types = "ccidccddd") |>
    rename(gene_id = phenotype_id) |>
    separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE,
             remove = FALSE) |>
    mutate(chrom = str_replace(chrom, "chr", "") |> as.integer(),
           pos = pos / 1e6)

tyr_tss <- tyr$pos[1] - tyr$tss_distance[1] / 1e6
tyr_thr <- perm |>
    filter(gene == "Tyr") |>
    pull(pval_nominal_threshold)

p1b <- tyr |>
    mutate(log10_pval = -log10(pval_nominal)) |>
    ggplot(aes(x = pos, y = log10_pval)) +
    geom_vline(xintercept = tyr_tss, color = "#777777") +
    geom_hline(yintercept = -log10(tyr_thr), color = "blue", lty = 2) +
    geom_point(size = 0.3) +
    geom_vline(xintercept = 151.709252, color = "#ff821c", lty = 3) +
    xlab("chr1 position (Mbp)") +
    ylab(expression(-log[10]*"p-value")) +
    theme_minimal() +
    annotate("text", x = 151.05, y = 3.2, color = "#666666", label = "Tyr\nTSS",
             hjust = 1) +
    annotate("text", x = 151.779252, y = 3.6, color = "#ff821c", label = "IOP QTL",
             hjust = 0) +
    annotate("text", x = 149, y = 3.7, color = "blue", label = "Tyr eQTL threshold",
             hjust = 0) +
    coord_cartesian(xlim = c(148.709252, 154.709252), expand = FALSE)
p1b

##########################
## Ndufaf6 associations ##
##########################

ndufaf6 <- read_tsv("inspect_qtl/ENSRNOG00000040040.4mbp.cis_qtl_pairs.txt",
                    col_types = "ccidccddd") |>
    rename(gene_id = phenotype_id) |>
    separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE,
             remove = FALSE) |>
    mutate(chrom = str_replace(chrom, "chr", "") |> as.integer(),
           pos = pos / 1e6)

ndufaf6_tss <- ndufaf6$pos[1] - ndufaf6$tss_distance[1] / 1e6
ndufaf6_thr <- perm |>
    filter(gene == "Ndufaf6") |>
    pull(pval_nominal_threshold)

p1c <- ndufaf6 |>
    mutate(log10_pval = -log10(pval_nominal)) |>
    ggplot(aes(x = pos, y = log10_pval)) +
    geom_vline(xintercept = ndufaf6_tss, color = "#777777") +
    geom_hline(yintercept = -log10(ndufaf6_thr), color = "blue", lty = 2) +
    geom_point(size = 0.3) +
    geom_vline(xintercept = 23.768259, color = "#e036e3", lty = 3) +
    xlab("chr5 position (Mbp)") +
    ylab(expression(-log[10]*"p-value")) +
    theme_minimal() +
    annotate("text", x = 24.4, y = 3, color = "#666666", label = "Ndufaf6\nTSS",
             hjust = 0) +
    annotate("text", x = 23.8, y = 3, color = "#e036e3", label = "IOP QTL",
             hjust = 0) +
    annotate("text", x = 21, y = 4.05, color = "blue", label = "Ndufaf6 eQTL threshold",
             hjust = 0) +
    coord_cartesian(xlim = c(20.768259, 26.768259), expand = FALSE)
p1c


###########################################
## Top assoc for genes near each IOP QTL ##
###########################################

# plot_qtl_region <- function(chr, posn) {
#     iop_color <- c("16" = "red", "1" = "#ff821c", "5" = "#e036e3")[as.character(chr)]
#     perm |>
#         mutate(pos = tss) |>
#         filter(chrom == chr,
#                pos > posn - 3.1,
#                pos < posn + 3.1) |>
#         mutate(pval = -log10(pval_nominal),
#                threshold = -log10(pval_nominal_threshold),
#                signif = pval_nominal < pval_nominal_threshold) |>
#         ggplot(aes(x = pos, y = pval, label = gene)) +
#         geom_segment(aes(x = pos, y = pval,
#                          xend = pos, yend = threshold),
#                      size = 0.3, alpha = 0.5, color = "blue") +
#         geom_vline(xintercept = posn, color = iop_color, lty = 3) +
#         geom_point(size = 1) +
#         annotate("text", x = posn + 0.07, y = 0.5, color = iop_color,
#                  label = "IOP QTL", hjust = 0) +
#         expand_limits(y = 0) +
#         coord_cartesian(xlim = c(posn - 3, posn + 3), expand = FALSE) +
#         xlab(str_glue("chr{chr} position (Mb)")) +
#         ylab(expression(-log[10]*"p-value")) +
#         theme_minimal()
# }

plot_qtl_region_labels <- function(chr, posn) {
    iop_color <- c("16" = "red", "1" = "#ff821c", "5" = "#e036e3")[as.character(chr)]
    perm |>
        mutate(pos = tss) |>
        filter(chrom == chr,
               pos > posn - 2,
               pos < posn + 2) |>
        mutate(pval = -log10(pval_nominal),
               threshold = -log10(pval_nominal_threshold),
               signif = pval_nominal < pval_nominal_threshold) |>
               # gene = if_else(pval_nominal < 0.001 | gene %in% c("Angpt2", "Tyr", "Ndufaf6"),
               #                gene, "")) |>
        ggplot(aes(x = pos, y = pval, label = gene)) +
        geom_segment(aes(x = pos, y = pval,
                         xend = pos, yend = threshold),
                     size = 0.3, alpha = 0.5) + #, color = "blue") +
        geom_vline(xintercept = posn, color = iop_color, lty = 3) +
        geom_point(size = 1) +
        # geom_point(aes(y = log10_threshold), shape = "-", size = 3, color = "red") +
        geom_segment(aes(x = pos - 0.075, xend = pos + 0.075, y = threshold,
                         yend = threshold), size = 0.25, color = "blue") +
        ggrepel::geom_text_repel(size = 3, segment.size = 0.2, #force_pull = 1e3,
                                 max.iter = 1e7, nudge_y = -0.3, min.segment.length = 0.2,
                                 max.overlaps = 100) +
        # geom_text(nudge_y = -0.5, size = 3) +
        annotate("text", x = posn + 0.05, y = 0.25, color = iop_color,
                 label = "IOP QTL", hjust = 0) +
        expand_limits(y = 0) +
        coord_cartesian(xlim = c(posn - 2, posn + 2), expand = FALSE) +
        xlab(str_glue("chr{chr} position (Mb)")) +
        ylab(expression(-log[10]*"p-value")) +
        theme_minimal()
}

# The labels are arranged badly. So make a version with and without them, and use
# the labeled version as a guide to annotate manually.

# p2 <- plot_qtl_region(16, 76.528844)
# p3 <- plot_qtl_region(1, 151.709252)
# p4 <- plot_qtl_region(5, 23.768259)
# 
# p1 / p2 / p3 / p4 + plot_layout(heights = c(1, 1, 1, 1))
# ggsave("figures/eQTL_figure.unlab.png", width = 8, height = 8)

p2 <- plot_qtl_region_labels(16, 76.528844)
p3 <- plot_qtl_region_labels(1, 151.709252)
p4 <- plot_qtl_region_labels(5, 23.768259)

p1 / p2 / p3 / p4 +
    plot_layout(heights = c(1, 1, 1, 1)) +
    plot_annotation(tag_levels = "A")
ggsave("figures/eQTL_figure.png", width = 8, height = 8)

# # Also make one with panels for Tyr and Ndufaf6:
# p1 / p2 / p1b / p3 / p1c / p4 + plot_layout(heights = c(1, 1, 1, 1, 1, 1))
# ggsave("figures/eQTL_figure.6panels.unlab.png", width = 8, height = 12)
