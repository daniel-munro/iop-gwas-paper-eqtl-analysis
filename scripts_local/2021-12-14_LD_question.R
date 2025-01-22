# Oksana asks: "What is the LD with IOP top SNP on chr 5 for Rbm12b and LOC679087?"

library(VariantAnnotation)
library(GenomicRanges)
library(tidyverse)

geno <- function(filename, ids) {
    snps <- tibble(variant_id = ids) |>
        separate(variant_id, c("chrom", "pos"), convert = TRUE) |>
        mutate(chrom = str_replace(chrom, "chr", ""))
    rng <- with(snps, GRanges(chrom, IRanges(pos, pos)))
    gt <- readGT(filename, param = rng)
    geno <- apply(gt, 2, function(x) c("0|0" = 0, "0|1" = 1, "1|0" = 1, "1|1" = 2)[x])
    rownames(geno) <- rownames(gt)
    geno
}

perm <- read_tsv("data/Eye.top_assoc.txt", col_types = "cc-c----------d--") |>
    filter(gene_name %in% c("Rbm12b", "LOC679087", "Ctsc"))

esnp <- perm |>
    select(gene_name, variant_id) |>
    deframe()

genos <- geno("data/genotype/eyes.vcf.gz",
              c(esnp[c("Rbm12b", "LOC679087")], "chr5:23768259"))
cor(t(genos)) ^ 2

# My own question, what about Ctsc?

genos2 <- geno("data/genotype/eyes.vcf.gz",
              c(esnp["Ctsc"], "chr1:151709252"))
cor(t(genos2)) ^ 2
