#!/usr/bin/env bash
# set -euo pipefail

# Genes of interest for IOP QTLs:
# GENE=ENSRNOG00000016696 # Angpt2
# GENE=ENSRNOG00000016421 # Tyr
# GENE=ENSRNOG00000040040 # Ndufaf6
GENE=ENSRNOG00000030719 # Csmd1
BED=data/tensorqtl/eyes.expression.bed.gz
NEWBED=data/tensorqtl/one_gene/$GENE.bed
PREFIX=$GENE
COVAR=data/tensorqtl/main4.PEER_covariates.txt
GENO=data/genotype/eyes
OUTDIR=data/tensorqtl/one_gene

# zcat $BED | head -1 > $NEWBED
# zcat $BED | grep $GENE >> $NEWBED
# bgzip $NEWBED

# Larger window for figure:
python -m tensorqtl --mode cis_nominal \
    $GENO \
    $NEWBED.gz \
    $PREFIX.4mbp \
    --covariates $COVAR \
    --window 4200000 \
    --output_dir $OUTDIR

