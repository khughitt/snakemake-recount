"""
Dimension reduction: PCA
kh (mar.2020)
"""
import pandas as pd
from sklearn.decomposition import PCA

# load gene expression data
dat = pd.read_feather(snakemake.input[0]).set_index('ensgene')

sample_ids = dat.columns

# convert to ndarray
dat = dat.to_numpy()

# perform sample-wise pca projections
# note that by default, PCA reduces column dimensionality, while KPCA/FastICA reduce row
# dimensionality
pca = PCA(n_components=snakemake.config['dimension_reduction']['pca']['num_dims'])
fit = pca.fit(dat)

res = pd.DataFrame(fit.components_)

res.columns = sample_ids
res.index = ['PC' + str(i) for i in range(1, res.shape[0] + 1)]

var_explained = pd.DataFrame(fit.explained_variance_ratio_)
var_explained.columns = ["ratio_var_explained"]
var_explained.index = res.index

# store projected data and variance explained
res.reset_index().rename(columns={"index": "dim"}).to_feather(snakemake.output[0])
var_explained.to_csv(snakemake.output[1])
