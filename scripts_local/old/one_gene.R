library(tidyverse)

cis <- read_tsv("inspect_qtl/ENSRNOG00000016696.cis_qtl_pairs.txt",
                col_types = "ccidiiddd") %>%
    arrange(pval_nominal)

write_tsv(cis, "inspect_qtl/Angpt2.cis_pairs.txt")
