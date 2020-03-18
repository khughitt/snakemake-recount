"""
Dimension reduction: PCA
KH (March 2020)
"""
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

# load gene expression data
dat = pd.read_csv(snakemake.input[0], sep='\t')

# drop gene id column convert to numpy array
dat = dat.iloc[:, :-1].to_numpy()

# perform sample-wise pca projections
pca = PCA(n_components=snakemake.config['dimension_reduction']['num_dims'])
fit = pca.fit(dat)

res = pd.DataFrame(fit.components_)

res.columns = dat.columns[:-1]
res.index = ['PC' + str(i) for i in range(1, res.shape[0] + 1)]

# store projected data and variance explained
np.savetxt(snakemake.output['fit'], res, delimiter='\t')
np.savetxt(snakemake.output['var'], fit.explained_variance_ratio_, delimiter='\t')

