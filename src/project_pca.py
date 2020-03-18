"""
Dimension reduction: PCA
KH (March 2020)
"""
import pandas as pd
from sklearn.decomposition import PCA

# load gene expression data
dat = pd.read_csv(snakemake.input[0], sep='\t')

sample_ids = dat.columns[:-1]

# drop gene id column convert to numpy array
dat = dat.iloc[:, :-1].to_numpy()

# perform sample-wise pca projections

# note that by default, PCA reduces column dimensionality, while KPCA/FastICA reduce row
# dimensionality
pca = PCA(n_components=snakemake.config['dimension_reduction']['num_dims'])
fit = pca.fit(dat)

res = pd.DataFrame(fit.components_)

res.columns = sample_ids
res.index = ['PC' + str(i) for i in range(1, res.shape[0] + 1)]

var_explained = pd.DataFrame(fit.explained_variance_ratio_)
var_explained.columns = ["ratio_var_explained"]
var_explained.index = res.index

# store projected data and variance explained
res.to_csv(snakemake.output['fit'], sep = '\t')
var_explained.to_csv(snakemake.output['var'], sep='\t')

