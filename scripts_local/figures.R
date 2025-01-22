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

perm <- read_tsv("data/Eye.cis_qtl.txt.gz", col_types = "ci----ci---d---ddd") |>
    rename(gene_id = phenotype_id) |>
    left_join(genes, by = "gene_id") |>
    # In case of duplicate genes, keep lowest p-value
    group_by(gene, chrom, tss) |>
    arrange(pval_nominal) |>
    slice(1) |>
    ungroup()

tss <- genes |>
    filter(gene %in% c("Tyr", "Ctsc", "Ndufaf6", "Plekhf2", "Csmd1", "Angpt2")) |>
    select(gene, tss) |>
    deframe()

thr <- perm |>
    filter(gene %in% c("Tyr", "Ctsc", "Ndufaf6", "Plekhf2", "Csmd1", "Angpt2")) |>
    mutate(thresh = -log10(pval_nominal_threshold)) |>
    select(gene, thresh) |>
    deframe()

pairs <- tibble(gene_id = c("ENSRNOG00000016696", "ENSRNOG00000016496", "ENSRNOG00000016421",
                            "ENSRNOG00000026662", "ENSRNOG00000040040", "ENSRNOG00000030719")) |>
    group_by(gene_id) |>
    summarise(
        read_tsv(str_glue("data/inspect_qtl/{gene_id}.4mbp.cis_qtl_pairs.txt.gz"),
                 col_types = "ccidccddd"),
        .groups = "drop"
    ) |>
    filter(!is.na(pval_nominal)) |>
    separate(variant_id, c("chrom", "pos"), sep = ":", convert = TRUE) |>
    mutate(pos = pos / 1e6,
           log10_pval = -log10(pval_nominal))

iop <- c("1" = 151.709252, "5" = 23.768259, "16" = 76.528844)
iop_colors <- c("1" = "#ff821c", "5" = "#e036e3", "16" = "red")
# window <- 1.5 # Mb on each side of IOP
xmin <- c("1" = 150.5, "5" = 23, "16" = 75.5)
xmax <- c("1" = 153.5, "5" = 26, "16" = 78)

######################
## Tyr associations ##
######################

tmp2c <- pairs |>
    filter(gene_id == "ENSRNOG00000016421")
ggplot(tmp2c, aes(x = pos, y = log10_pval)) +
    geom_vline(xintercept = tss["Tyr"], color = "#777777") +
    # geom_hline(yintercept = thr["Tyr"], color = "blue", lty = 2) +
    geom_segment(aes(x = x, xend = xend, y = y, yend = yend),
                 data = tibble(x = tss["Tyr"] - 0.5, xend = tss["Tyr"] + 0.5,
                               y = thr["Tyr"], yend = thr["Tyr"]),
                 color = "blue", lty = 2) +
    geom_point(size = 0.3) +
    geom_vline(xintercept = iop["1"], color = iop_colors["1"], lty = 3) +
    xlab("chr1 position (Mb)") +
    ylab(expression(-log[10](P))) +
    geom_richtext(data = tibble(a=1), x = tss["Tyr"] + 0.0, y = 2.75, color = "#666666",
                  label = "*Tyr* TSS", fill = NA, label.color = NA, hjust = 1) +
    annotate("text", x = iop["1"] + 0.03, y = 2.75, color = iop_colors["1"],
             label = "IOP QTL", hjust = 0) +
    # geom_richtext(data = tibble(a=1), x = 152.5, y = 3.55, color = "blue",
    #               label = "*Tyr* eQTL threshold", fill = NA, label.color = NA, hjust = 0) +
    geom_richtext(data = tibble(a=1), x = tss["Tyr"] + 0.02, y = 3.55, color = "blue",
                  label = "cis-eQTL threshold", fill = NA, label.color = NA, hjust = 0) +
    coord_cartesian(xlim = c(xmin["1"], xmax["1"]), expand = FALSE) +
    expand_limits(y = max(tmp2c$log10_pval) + 0.2) +
    theme_bw() +
    theme(plot.margin = unit(c(5.5, 20, 5.5, 5.5), "pt"),
          panel.grid = element_blank())
# ggsave("figures/eQTL_figure.2c.png", width = 8, height = 2, bg = "white")

#######################
## Ctsc associations ##
#######################

tmp2c <- pairs |>
    filter(gene_id == "ENSRNOG00000016496")
ggplot(tmp2c, aes(x = pos, y = log10_pval)) +
    geom_vline(xintercept = tss["Ctsc"], color = "#777777") +
    geom_segment(aes(x = x, xend = xend, y = y, yend = yend),
                 data = tibble(x = tss["Ctsc"] - 0.5, xend = tss["Ctsc"] + 0.5,
                               y = thr["Ctsc"], yend = thr["Ctsc"]),
                 color = "blue", lty = 2) +
    geom_point(size = 0.3) +
    geom_vline(xintercept = iop["1"], color = iop_colors["1"], lty = 3) +
    xlab("chr1 position (Mb)") +
    ylab(expression(-log[10](P))) +
    geom_richtext(data = tibble(a=1), x = tss["Ctsc"] + 0.0, y = 2.2, color = "#666666",
                  label = "*Ctsc* TSS", fill = NA, label.color = NA, hjust = 0) +
    annotate("text", x = iop["1"] - 0.03, y = 2.2, color = iop_colors["1"],
             label = "IOP QTL", hjust = 1) +
    geom_richtext(data = tibble(a=1), x = tss["Ctsc"] + 0.53, y = 3.37, color = "blue",
                  label = "cis-eQTL threshold", fill = NA, label.color = NA, hjust = 0) +
    coord_cartesian(xlim = c(xmin["1"], xmax["1"]), expand = FALSE) +
    expand_limits(y = max(tmp2c$log10_pval) + 0.2) +
    theme_bw() +
    theme(plot.margin = unit(c(5.5, 20, 5.5, 5.5), "pt"),
          panel.grid = element_blank())
ggsave("figures/eQTL_figure.2c.png", width = 8, height = 2, bg = "white")

##########################
## Ndufaf6 associations ##
##########################

tmp3c <- pairs |>
    filter(gene_id == "ENSRNOG00000040040")
ggplot(tmp3c, aes(x = pos, y = log10_pval)) +
    geom_vline(xintercept = tss["Ndufaf6"], color = "#777777") +
    # geom_hline(yintercept = thr["Ndufaf6"], color = "blue", lty = 2) +
    geom_segment(aes(x = x, xend = xend, y = y, yend = yend),
                 data = tibble(x = tss["Ndufaf6"] - 0.5, xend = tss["Ndufaf6"] + 0.5,
                               y = thr["Ndufaf6"], yend = thr["Ndufaf6"]),
                 color = "blue", lty = 2) +
    geom_point(size = 0.3) +
    geom_vline(xintercept = iop["5"], color = iop_colors["5"], lty = 3) +
    xlab("chr5 position (Mb)") +
    ylab(expression(-log[10](P))) +
    geom_richtext(data = tibble(a=1), x = tss["Ndufaf6"] + 0.0, y = 3, color = "#666666",
                  label = "*Ndufaf6* TSS", fill = NA, label.color = NA, hjust = 0) +
    annotate("text", x = iop["5"] + 0.03, y = 3, color = iop_colors["5"],
             label = "IOP QTL", hjust = 0) +
    # geom_richtext(data = tibble(a=1), x = 25, y = 3.25, color = "blue",
    #               label = "*Ndufaf6* eQTL threshold", fill = NA, label.color = NA, hjust = 0) +
    geom_richtext(data = tibble(a=1), x = 24.85, y = 3.5, color = "blue",
                  label = "cis-eQTL threshold", fill = NA, label.color = NA, hjust = 0) +
    coord_cartesian(xlim = c(xmin["5"], xmax["5"]), expand = FALSE) +
    expand_limits(y = max(c(tmp3c$log10_pval, thr["Ndufaf6"])) + 0.2) +
    theme_bw() +
    theme(plot.margin = unit(c(5.5, 20, 5.5, 5.5), "pt"),
          panel.grid = element_blank())
ggsave("figures/eQTL_figure.3c.png", width = 8, height = 2, bg = "white")

##########################
## Plekhf2 associations ##
##########################

tmpPlekhf2 <- pairs |>
    filter(gene_id == "ENSRNOG00000026662",
           pos > 23 - 0.1, pos < 26 + 0.1)
ggplot(tmpPlekhf2, aes(x = pos, y = log10_pval)) +
    geom_vline(xintercept = tss["Plekhf2"], color = "#777777") +
    # geom_hline(yintercept = thr["Plekhf2"], color = "blue", lty = 2) +
    geom_segment(aes(x = x, xend = xend, y = y, yend = yend),
                 data = tibble(x = tss["Plekhf2"] - 0.5, xend = tss["Plekhf2"] + 0.5,
                               y = thr["Plekhf2"], yend = thr["Plekhf2"]),
                 color = "blue", lty = 2) +
    geom_point(size = 0.3) +
    geom_vline(xintercept = iop["5"], color = iop_colors["5"], lty = 3) +
    xlab("chr5 position (Mb)") +
    ylab(expression(-log[10](P))) +
    geom_richtext(data = tibble(a=1), x = tss["Plekhf2"] + 0.0, y = 0.3, color = "#666666",
                  label = "*Plekhf2* TSS", fill = NA, label.color = NA, hjust = 0) +
    annotate("text", x = iop["5"] + 0.03, y = 3, color = iop_colors["5"],
             label = "IOP QTL", hjust = 0) +
    # geom_richtext(data = tibble(a=1), x = 25, y = 3.25, color = "blue",
    #               label = "*Plekhf2* eQTL threshold", fill = NA, label.color = NA, hjust = 0) +
    geom_richtext(data = tibble(a=1), x = 24.78, y = 3.5, color = "blue",
                  label = "cis-eQTL threshold", fill = NA, label.color = NA, hjust = 0) +
    coord_cartesian(xlim = c(xmin["5"], xmax["5"]), expand = FALSE) +
    expand_limits(y = max(c(tmpPlekhf2$log10_pval, thr["Plekhf2"])) + 0.3) +
    theme_bw() +
    theme(plot.margin = unit(c(5.5, 20, 5.5, 5.5), "pt"),
          panel.grid = element_blank())
ggsave("figures/eQTL_Plekhf2.png", width = 8, height = 2, bg = "white")

#########################
## Csmd1 associations ##
#########################

pairs |>
    filter(gene_id == "ENSRNOG00000030719") |>
    ggplot(aes(x = pos, y = log10_pval)) +
    geom_vline(xintercept = tss["Csmd1"], color = "#777777") +
    # geom_hline(yintercept = thr["Csmd1"], color = "blue", lty = 2) +
    geom_segment(aes(x = x, xend = xend, y = y, yend = yend),
                 data = tibble(x = tss["Csmd1"] - 0.5, xend = tss["Csmd1"] + 0.5,
                               y = thr["Csmd1"], yend = thr["Csmd1"]),
                 color = "blue", lty = 2) +
    geom_point(size = 0.3) +
    geom_vline(xintercept = iop["16"], color = iop_colors["16"], lty = 3) +
    xlab("chr16 position (Mb)") +
    ylab(expression(-log[10](P))) +
    geom_richtext(data = tibble(a=1), x = tss["Csmd1"] + 0.0, y = 3, color = "#666666",
                  label = "*Csmd1* TSS", fill = NA, label.color = NA, hjust = 0) +
    annotate("text", x = iop["16"] + 0.03, y = 2.5, color = iop_colors["16"],
             label = "IOP QTL", hjust = 0) +
    # geom_richtext(data = tibble(a=1), x = xmin["16"] + 0.2, y = 3.1, color = "blue",
    #               label = "*Csmd1* eQTL threshold", fill = NA, label.color = NA, hjust = 0) +
    geom_richtext(data = tibble(a=1), x = 77.25, y = thr["Csmd1"] + 0.2, color = "blue",
                  label = "cis-eQTL threshold", fill = NA, label.color = NA, hjust = 0) +
    coord_cartesian(xlim = c(xmin["16"], xmax["16"]), expand = FALSE) +
    expand_limits(y = 4.2) +
    theme_bw() +
    theme(plot.margin = unit(c(5.5, 20, 5.5, 5.5), "pt"),
          panel.grid = element_blank())
ggsave("figures/eQTL_figure.4c.png", width = 8, height = 2, bg = "white")

#########################
## Angpt2 associations ##
#########################

tmpAngpt2 <- pairs |>
    filter(gene_id == "ENSRNOG00000016696")
ggplot(tmpAngpt2, aes(x = pos, y = log10_pval)) +
    geom_vline(xintercept = tss["Angpt2"], color = "#777777") +
    geom_hline(yintercept = thr["Angpt2"], color = "blue", lty = 2) +
    geom_point(size = 0.3) +
    geom_vline(xintercept = iop["16"], color = iop_colors["16"], lty = 3) +
    xlab("chr16 position (Mb)") +
    ylab(expression(-log[10](P))) +
    geom_richtext(data = tibble(a=1), x = tss["Angpt2"] + 0.0, y = 3, color = "#666666",
                  label = "*Angpt2* TSS", fill = NA, label.color = NA, hjust = 0) +
    annotate("text", x = iop["16"] + 0.03, y = 3.5, color = iop_colors["16"],
             label = "IOP QTL", hjust = 0) +
    geom_richtext(data = tibble(a=1), x = 76.5, y = thr["Angpt2"] + 0.2, color = "blue",
                  label = "*Angpt2* eQTL threshold", fill = NA, label.color = NA, hjust = 0) +
    coord_cartesian(xlim = c(xmin["16"] - 2, xmax["16"]), expand = FALSE) +
    expand_limits(y = max(c(tmpAngpt2$log10_pval, thr["Angpt2"])) + 0.2) +
    theme_bw() +
    theme(plot.margin = unit(c(5.5, 20, 5.5, 5.5), "pt"),
          panel.grid = element_blank())
ggsave("figures/eQTL_Angpt2.png", width = 8, height = 2, bg = "white")


###########################################
## Top assoc for genes near each IOP QTL ##
###########################################

plot_qtl_region_labels <- function(chr) {
    tmp <- perm |>
        filter(chrom == chr,
               # tss > posn - window,
               # tss < posn + window) |>
               tss > xmin[chr],
               tss < xmax[chr]) |>
        mutate(pval = -log10(pval_nominal),
               threshold = -log10(pval_nominal_threshold),
               signif = pval_nominal < pval_nominal_threshold)
    IOP_label_y <- c("1" = 5, "5" = 5, "16" = 3)[chr]
    ggplot(tmp, aes(x = tss, y = pval, label = gene)) +
        geom_segment(aes(x = tss, y = pval, xend = tss, yend = threshold),
                     size = 0.2, color = "#555555") +
        # geom_point(aes(y = log10_threshold), shape = "-", size = 3, color = "red") +
        geom_segment(aes(x = tss - 0.05, xend = tss + 0.05, y = threshold,
                         yend = threshold), size = 0.3, color = "blue") +
        geom_point(size = 1) +
        geom_vline(xintercept = iop[chr], color = iop_colors[chr], lty = 3) +
        ggrepel::geom_text_repel(size = 3, segment.size = 0.2, #force_pull = 1e3,
                                 max.iter = 1e7, nudge_y = -0.3, min.segment.length = 0,
                                 max.overlaps = 100, fontface = "italic") +
        # geom_text(nudge_y = -0.5, size = 3) +
        annotate("text", x = iop[chr] + 0.03, y = IOP_label_y, color = iop_colors[chr],
                 label = "IOP QTL", hjust = 0) +
        # coord_cartesian(xlim = c(posn - window, posn + window), expand = FALSE) +
        coord_cartesian(xlim = c(xmin[chr], xmax[chr]), expand = FALSE) +
        expand_limits(y = c(0, max(tmp$pval) + 0.2)) +
        xlab(str_glue("chr{chr} position (Mb)")) +
        ylab(expression(-log[10](P))) +
        theme_bw() +
        theme(plot.margin = unit(c(5.5, 20, 5.5, 5.5), "pt"),
              panel.grid = element_blank())
}

plot_qtl_region_labels("1")
ggsave("figures/eQTL_figure.2b.svg", width = 8, height = 2, bg = "white")
ggsave("figures/eQTL_figure.2b.guide.png", width = 8, height = 2, bg = "white")

plot_qtl_region_labels("5")
ggsave("figures/eQTL_figure.3b.svg", width = 8, height = 2, bg = "white")
ggsave("figures/eQTL_figure.3b.guide.png", width = 8, height = 2, bg = "white")

plot_qtl_region_labels("16")
ggsave("figures/eQTL_figure.4b.svg", width = 8, height = 2, bg = "white")
ggsave("figures/eQTL_figure.4b.guide.png", width = 8, height = 2, bg = "white")

##################
## Effect plots ##
##################

# Input these pairs into https://ratgtex.org/eqtl-dashboard/
perm |>
    filter(gene_id %in% pairs$gene_id) |>
    mutate(pair = str_c(gene_id, variant_id, sep = ",")) |>
    select(pair)
