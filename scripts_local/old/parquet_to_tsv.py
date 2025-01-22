import sys
import pandas as pd

# pd.read_parquet("chr16_76543372.cis_qtl_pairs.16.parquet").to_csv("chr16_76543372.cis_qtl_pairs.txt", sep="\t", index=False)

# pd.read_parquet("chr16_76543372.trans_qtl_pairs.parquet").to_csv("chr16_76543372.trans_qtl_pairs.txt", sep="\t", index=False)

# pd.read_parquet("ENSRNOG00000002227.cis_qtl_pairs.14.parquet").to_csv("ENSRNOG00000002227.cis_qtl_pairs.txt", sep="\t", index=False)

pd.read_parquet(sys.argv[1]).to_csv(sys.argv[2], sep="\t", index=False)
