library(tidyverse)
library(patchwork)
library(ggtext)

genes <- read_tsv("data/reference/GENES_RAT.txt",
                  col_types = cols(ENSEMBL_ID = "c", SYMBOL = "c", CHROMOSOME_6.0 = "c",
                                   START_POS_6.0 = "c", STOP_POS_6.0 = "c",
                                   STRAND_6.0 = "c", .default = "-"),
                  skip = 83) |>
    rename(gene_id = ENSEMBL_ID,
           gene = SYMBOL,
           chrom = CHROMOSOME_6.0,
           start = START_POS_6.0,
           stop = STOP_POS_6.0,
           strand = STRAND_6.0) |>
    separate_rows(gene_id, sep = ";") |>
    filter(!is.na(gene_id),
           !is.na(start)) |>
    # filter(!str_detect(gene, "^LOC\\d+$")) |>
    group_by(gene_id) |>
    slice(1) |>
    ungroup() |>
    # For genes with multiple locations, keep first:
    mutate(chrom = str_split(chrom, ";", simplify = TRUE)[, 1],
           strand = str_split(strand, ";", simplify = TRUE)[, 1],
           start = str_split(start, ";", simplify = TRUE)[, 1] |> as.integer(),
           stop = str_split(stop, ";", simplify = TRUE)[, 1] |> as.integer(),
           tss = if_else(strand == "+", start, stop) / 1e6)

perm <- read_tsv("data/tensorqtl/main5.cis_qtl.txt.gz",
                 col_types = "ci----cd----d---ddd") |>
    rename(gene_id = phenotype_id) |>
    left_join(genes, by = "gene_id") |>
    # In case of duplicate genes, keep lowest p-value
    group_by(gene, chrom, tss) |>
    arrange(pval_nominal) |>
    slice(1) |>
    ungroup()

tss <- genes |>
    filter(gene %in% c("Angpt2", "Csmd1", "Tyr", "Ndufaf6")) |>
    select(gene, tss) |>
    deframe()

thr <- perm |>
    filter(gene %in% c("Angpt2", "Csmd1", "Tyr", "Ndufaf6")) |>
    mutate(thresh = -log10(pval_nominal_threshold)) |>
    select(gene, thresh) |>
    deframe()

pairs <- tibble(gene_id = c("ENSRNOG00000016696", "ENSRNOG00000030719",
                            "ENSRNOG00000016421", "ENSRNOG00000040040")) |>
    group_by(gene_id) |>
    summarise(
        read_tsv(str_glue("inspect_qtl/{gene_id}.4mbp.cis_qtl_pairs.txt"),
                 col_types = "ccidccddd"),
        .groups = "drop"
    ) |>
    filter(!is.na(pval_nominal)) |>
    separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE) |>
    mutate(pos = pos / 1e6,
           log10_pval = -log10(pval_nominal))

# Since different TSS was used for Csmd1, lowest p-vals in cis-window were missed:
perm <- perm |>
    mutate(pval_nominal = if_else(
        gene == "Csmd1",
        min(pairs$pval_nominal[pairs$gene_id == "ENSRNOG00000030719"]),
        pval_nominal
    ))

iop <- c("16" = 76.528844, "1" = 151.709252, "5" = 23.768259)
iop_colors <- c("16" = "red", "1" = "#ff821c", "5" = "#e036e3")
# window <- 1.5 # Mb on each side of IOP
xmin <- c("16" = 75.5, "1" = 150.5, "5" = 23)
xmax <- c("16" = 78, "1" = 153.5, "5" = 26)

# #########################
# ## Angpt2 associations ##
# #########################
# 
# p1 <- pairs |>
#     filter(gene_id == "ENSRNOG00000016696") |>
#     ggplot(aes(x = pos, y = log10_pval)) +
#     geom_vline(xintercept = tss["Angpt2"], color = "#777777") +
#     geom_hline(yintercept = thr["Angpt2"], color = "blue", lty = 2) +
#     geom_point(size = 0.3) +
#     geom_vline(xintercept = iop["16"], color = iop_colors["16"], lty = 3) +
#     xlab("chr16 position (Mb)") +
#     ylab(expression(-log[10]*"p-value")) +
#     theme_minimal() +
#     annotate("text", x = tss["Angpt2"] + 0.03, y = 3.25, color = "#666666",
#              label = "Angpt2 TSS", hjust = 0) +
#     annotate("text", x = iop["16"] + 0.03, y = 3, color = iop_colors["16"],
#              label = "IOP QTL", hjust = 0) +
#     annotate("text", x = iop["16"] - window + 0.2, y = 3.95, color = "blue",
#              label = "Angpt2 eQTL threshold", hjust = 0) +
#     coord_cartesian(xlim = c(iop["16"] - window, iop["16"] + window), expand = FALSE)
# p1
# 
# # ggsave("figures/Angpt2.png")

#########################
## Csmd1 associations ##
#########################

p1a <- pairs |>
    filter(gene_id == "ENSRNOG00000030719") |>
    ggplot(aes(x = pos, y = log10_pval)) +
    geom_vline(xintercept = tss["Csmd1"], color = "#777777") +
    geom_hline(yintercept = thr["Csmd1"], color = "blue", lty = 2) +
    geom_point(size = 0.3) +
    geom_vline(xintercept = iop["16"], color = iop_colors["16"], lty = 3) +
    xlab("chr16 position (Mb)") +
    ylab(expression(-log[10]*"p-value")) +
    # annotate("text", x = tss["Csmd1"] + 0.03, y = 3.25, color = "#666666",
    #          label = "Csmd1 TSS", hjust = 0) +
    geom_richtext(data = tibble(a=1), x = tss["Csmd1"] + 0.0, y = 3.25, color = "#666666",
                  label = "*Csmd1* TSS", fill = NA, label.color = NA, hjust = 0) +
    annotate("text", x = iop["16"] + 0.03, y = 2, color = iop_colors["16"],
             label = "IOP QTL", hjust = 0) +
    # annotate("text", x = iop["16"] - window + 0.2, y = 3.65, color = "blue",
    #          label = "Csmd1 eQTL threshold", hjust = 0) +
    geom_richtext(data = tibble(a=1), x = xmin["16"] + 0.2, y = 3.65, color = "blue",
              label = "*Csmd1* eQTL threshold", fill = NA, label.color = NA, hjust = 0) +
    # coord_cartesian(xlim = c(iop["16"] - window, iop["16"] + window), expand = FALSE) +
    coord_cartesian(xlim = c(xmin["16"], xmax["16"]), expand = FALSE) +
    expand_limits(y = 4.2) +
    theme_minimal() +
    theme(plot.margin = unit(c(5.5, 20, 5.5, 5.5), "pt"))
p1a

######################
## Tyr associations ##
######################

p1b <- pairs |>
    filter(gene_id == "ENSRNOG00000016421") |>
    ggplot(aes(x = pos, y = log10_pval)) +
    geom_vline(xintercept = tss["Tyr"], color = "#777777") +
    geom_hline(yintercept = thr["Tyr"], color = "blue", lty = 2) +
    geom_point(size = 0.3) +
    geom_vline(xintercept = iop["1"], color = iop_colors["1"], lty = 3) +
    xlab("chr1 position (Mb)") +
    ylab(expression(-log[10]*"p-value")) +
    # annotate("text", x = tss["Tyr"] - 0.07, y = 3.2, color = "#666666", label = "Tyr TSS",
    #          hjust = 1) +
    geom_richtext(data = tibble(a=1), x = tss["Tyr"] + 0.0, y = 3.25, color = "#666666",
                  label = "*Tyr* TSS", fill = NA, label.color = NA, hjust = 0) +
    annotate("text", x = iop["1"] + 0.03, y = 3.6, color = iop_colors["1"], label = "IOP QTL",
             hjust = 0) +
    # annotate("text", x = 149, y = 3.7, color = "blue", label = "Tyr eQTL threshold",
    #          hjust = 0) +
    geom_richtext(data = tibble(a=1), x = xmin["1"] + 0.02, y = 3.65, color = "blue",
                  label = "*Tyr* eQTL threshold", fill = NA, label.color = NA, hjust = 0) +
    # coord_cartesian(xlim = c(iop["1"] - window, iop["1"] + window), expand = FALSE)
    coord_cartesian(xlim = c(xmin["1"], xmax["1"]), expand = FALSE) +
    theme_minimal() +
    theme(plot.margin = unit(c(5.5, 20, 5.5, 5.5), "pt"))
p1b

##########################
## Ndufaf6 associations ##
##########################

p1c <- pairs |>
    filter(gene_id == "ENSRNOG00000040040") |>
    ggplot(aes(x = pos, y = log10_pval)) +
    geom_vline(xintercept = tss["Ndufaf6"], color = "#777777") +
    geom_hline(yintercept = thr["Ndufaf6"], color = "blue", lty = 2) +
    geom_point(size = 0.3) +
    geom_vline(xintercept = iop["5"], color = iop_colors["5"], lty = 3) +
    xlab("chr5 position (Mb)") +
    ylab(expression(-log[10]*"p-value")) +
    # annotate("text", x = tss["Ndufaf6"] + 0.07, y = 3, color = "#666666", label = "Ndufaf6 TSS",
    #          hjust = 0) +
    geom_richtext(data = tibble(a=1), x = tss["Ndufaf6"] + 0.0, y = 3, color = "#666666",
                  label = "*Ndufaf6* TSS", fill = NA, label.color = NA, hjust = 0) +
    annotate("text", x = iop["5"] + 0.03, y = 3, color = iop_colors["5"], label = "IOP QTL",
             hjust = 0) +
    # annotate("text", x = 21, y = 4.05, color = "blue", label = "Ndufaf6 eQTL threshold",
    #          hjust = 0) +
    geom_richtext(data = tibble(a=1), x = xmin["5"] + 0.05, y = 4, color = "blue",
                  label = "*Ndufaf6* eQTL threshold", fill = NA, label.color = NA, hjust = 0) +
    # coord_cartesian(xlim = c(iop["5"] - window, iop["5"] + window), expand = FALSE)
    coord_cartesian(xlim = c(xmin["5"], xmax["5"]), expand = FALSE) +
    theme_minimal() +
    theme(plot.margin = unit(c(5.5, 20, 5.5, 5.5), "pt"))
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

plot_qtl_region_labels <- function(chr) {
    perm |>
        filter(chrom == chr,
               # tss > posn - window,
               # tss < posn + window) |>
               tss > xmin[chr],
               tss < xmax[chr]) |>
        mutate(pval = -log10(pval_nominal),
               threshold = -log10(pval_nominal_threshold),
               signif = pval_nominal < pval_nominal_threshold) |>
        ggplot(aes(x = tss, y = pval, label = gene)) +
        geom_segment(aes(x = tss, y = pval, xend = tss, yend = threshold),
                     size = 0.2, color = "#555555") +
        # geom_point(aes(y = log10_threshold), shape = "-", size = 3, color = "red") +
        geom_segment(aes(x = tss - 0.05, xend = tss + 0.05, y = threshold,
                         yend = threshold), size = 0.3, color = "blue") +
        geom_point(size = 1) +
        ggrepel::geom_text_repel(size = 3, segment.size = 0.2, #force_pull = 1e3,
                                 max.iter = 1e7, nudge_y = -0.3, min.segment.length = 0.2,
                                 max.overlaps = 100, fontface = "italic") +
        # geom_text(nudge_y = -0.5, size = 3) +
        geom_vline(xintercept = iop[chr], color = iop_colors[chr], lty = 3) +
        annotate("text", x = iop[chr] + 0.03, y = 2, color = iop_colors[chr],
                 label = "IOP QTL", hjust = 0) +
        expand_limits(y = 0) +
        # coord_cartesian(xlim = c(posn - window, posn + window), expand = FALSE) +
        coord_cartesian(xlim = c(xmin[chr], xmax[chr]), expand = FALSE) +
        xlab(str_glue("chr{chr} position (Mb)")) +
        ylab(expression(-log[10]*"p-value")) +
        theme_minimal() +
        theme(plot.margin = unit(c(5.5, 20, 5.5, 5.5), "pt"))
}

# The labels are arranged badly. So make a version with and without them, and use
# the labeled version as a guide to annotate manually.

# p2 <- plot_qtl_region(16, 76.528844)
# p3 <- plot_qtl_region(1, 151.709252)
# p4 <- plot_qtl_region(5, 23.768259)
# 
# p1 / p2 / p3 / p4 + plot_layout(heights = c(1, 1, 1, 1))
# ggsave("figures/eQTL_figure.unlab.png", width = 8, height = 8)

p2a <- plot_qtl_region_labels("16")
p2b <- plot_qtl_region_labels("1")
p2c <- plot_qtl_region_labels("5")

# p1a / p2a / p2b / p2c +
#     plot_annotation(tag_levels = "A")
# ggsave("figures/eQTL_figure.png", width = 8, height = 8)

# # Also make one with panels for Tyr and Ndufaf6:
# p1 / p2 / p1b / p3 / p1c / p4 + plot_layout(heights = c(1, 1, 1, 1, 1, 1))
# ggsave("figures/eQTL_figure.6panels.unlab.png", width = 8, height = 12)

ggsave("figures/eQTL_figure.1.png", p1a, width = 8, height = 2, bg = "white")
ggsave("figures/eQTL_figure.2.png", p1b, width = 8, height = 2, bg = "white")
ggsave("figures/eQTL_figure.3.png", p1c, width = 8, height = 2, bg = "white")
ggsave("figures/eQTL_figure.4.png", p2a, width = 8, height = 2, bg = "white")
ggsave("figures/eQTL_figure.5.png", p2b, width = 8, height = 2, bg = "white")
ggsave("figures/eQTL_figure.6.png", p2c, width = 8, height = 2, bg = "white")
