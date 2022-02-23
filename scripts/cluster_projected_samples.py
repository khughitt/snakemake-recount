"""
Projected data sample clustering

Performs binary clustering independently along each dimension of the projected data
"""
import pandas as pd
from sklearn.cluster import KMeans

# load projected data
dat = pd.read_feather(snakemake.input[0]).set_index('dim')

# initialize KMeans engine
kmeans = KMeans(n_clusters=2, random_state=snakemake.config['random_seed'])

rows = []

# cluster along each dimension of projected data independently
num_dims = dat.shape[0]

for i in range(num_dims):
    mat = dat.iloc[i, :].to_numpy().reshape(-1, 1)
    rows.append(list(kmeans.fit(mat).labels_))

# create a dataframe of cluster assignments
res = pd.DataFrame(rows, index=dat.index, columns=dat.columns)

# store result
res.reset_index().rename(columns={"index": "dim"}).to_feather(snakemake.output[0])
