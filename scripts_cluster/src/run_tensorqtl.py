# Run in python script to use random_tiebreak option
# Based on examples at https://github.com/broadinstitute/tensorqtl

import argparse
import os
import pandas as pd
import tensorqtl
from tensorqtl import genotypeio, cis

parser = argparse.ArgumentParser(description="Run tensorQTL from command line but with extra options")
parser.add_argument("geno_prefix")
parser.add_argument("expression")
parser.add_argument("output", help="Output file")
parser.add_argument("--covariates", help="Covariates file")
parser.add_argument("--window", help="Window extent (on either side of TSS)", type=int, default=1000000)
parser.add_argument("--cis_output", help="For cis_independent mode, output of tensorQTL in cis mode")
parser.add_argument("--mode", required=True, help="Run mode: currently either cis or cis_independent")
args = parser.parse_args()

pheno, pheno_pos = tensorqtl.read_phenotype_bed(args.expression)
if args.covariates is not None:
    covar = pd.read_csv(args.covariates, sep="\t", index_col=0).T
else:
    covar = None
pr = genotypeio.PlinkReader(args.geno_prefix)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

if args.mode == "cis":
    d = cis.map_cis(genotype_df, variant_df, pheno, pheno_pos, covar, window=args.window, random_tiebreak=True)
    tensorqtl.calculate_qvalues(d, qvalue_lambda=0.85)
    d.to_csv(args.output, sep="\t", float_format="%.6g")
elif args.mode == "cis_independent":
    cis_df = pd.read_csv(args.cis_output, sep="\t", index_col=0)
    d = cis.map_independent(genotype_df, variant_df, cis_df, pheno, pheno_pos, covar, window=args.window, random_tiebreak=True)
    d.to_csv(args.output, sep="\t", index=False, float_format="%.6g")
else:
    print("Mode not recognized.")
