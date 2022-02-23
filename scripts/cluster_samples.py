"""
Filtered data sample clustering

Performs binary using non-projected version of data
"""
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

# load data
dat = pd.read_feather(snakemake.input[0]).set_index('ensgene')
sample_ids = dat.columns

# scale data
scaler = StandardScaler()
dat = scaler.fit_transform(dat)

# cluster samples
kmeans = KMeans(n_clusters=2, random_state=snakemake.config['random_seed'])
cluster_labels = list(kmeans.fit(dat.T).labels_)

# store result
res = pd.DataFrame({'sample_id': sample_ids, 'kmeans': cluster_labels})
res.to_feather(snakemake.output[0])
