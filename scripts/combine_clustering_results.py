"""
Combines clustering results into a single dataframe
"""
import pandas as pd

# clusters for original data are stored as a <sample x 2> dataframe
combined = pd.read_feather(snakemake.input[0]).set_index('sample_id')

# clusters for projected data are stored as <n x sample> dataframes
for infile in snakemake.input[1:]:
    dat = pd.read_feather(infile).set_index('dim').T
    combined = pd.concat([combined, dat], axis=1)

combined.reset_index().rename(columns={"index": "sample_id"}).to_feather(snakemake.output[0])
