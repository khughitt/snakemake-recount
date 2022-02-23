"""
Dimension reduction: PCA
kh (mar.2020)
"""
import pandas as pd
from sklearn.decomposition import KernelPCA

# load gene expression data
dat = pd.read_feather(snakemake.input[0]).set_index('ensgene')

sample_ids = dat.columns

# drop gene id column convert to numpy array
dat = dat.to_numpy()

# perform sample-wise PCA projections;
# note that by default, PCA reduces column dimensionality, while KPCA/FastICA reduce row
# dimensionality
ind = "kpca_" + snakemake.params['kernel']
num_dims = snakemake.config['dimension_reduction'][ind]['num_dims']

transformer = KernelPCA(n_components=num_dims, kernel=snakemake.params['kernel'],
                        random_state=snakemake.config['random_seed'])

fit = transformer.fit_transform(dat.T)

# convert result to a dataframe with PCs as row and samples as columns
res = pd.DataFrame(fit.T)

res.columns = sample_ids
res.index = ['KPC' + str(i) for i in range(1, res.shape[0] + 1)]

# store projected data and variance explained
res.reset_index().rename(columns={"index": "dim"}).to_feather(snakemake.output[0])
