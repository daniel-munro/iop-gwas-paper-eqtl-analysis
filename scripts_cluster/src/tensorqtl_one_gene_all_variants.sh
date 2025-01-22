#!/usr/bin/env bash
# set -euo pipefail

## Genes of interest for IOP QTLs:
# GENE=ENSRNOG00000016696 # Angpt2
# GENE=ENSRNOG00000016421 # Tyr
# GENE=ENSRNOG00000040040 # Ndufaf6
# GENE=ENSRNOG00000030719 # Csmd1
# GENE=ENSRNOG00000016496 # Ctsc
GENE=ENSRNOG00000026662 # Plekhf2
BED=data/Eye.expr.iqn.filtered.bed.gz
NEWBED=data/inspect_qtl/$GENE.bed
PREFIX=$GENE
COVAR=data/covar.txt
GENO=data/genotype/geno
OUTDIR=data/inspect_qtl

echo $GENE

# zcat $BED | head -1 > $NEWBED
# zcat $BED | grep $GENE >> $NEWBED
# bgzip $NEWBED

# # Larger window for figure:
# python -m tensorqtl --mode cis_nominal \
#     $GENO \
#     $NEWBED.gz \
#     $PREFIX.4mbp \
#     --covariates $COVAR \
#     --window 4000000 \
#     --output_dir $OUTDIR

PARQUET=$OUTDIR/$PREFIX.4mbp.cis_qtl_pairs.*.parquet
PARQUET2=$OUTDIR/$PREFIX.4mbp.cis_qtl_pairs.parquet
mv -i $PARQUET $PARQUET2
OUTFILE=$OUTDIR/$PREFIX.4mbp.cis_qtl_pairs.txt.gz
python3 -c "import pandas as pd; pd.read_parquet('$PARQUET2').to_csv('$OUTFILE', sep='\t', index=False)"
