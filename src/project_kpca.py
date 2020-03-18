"""
Dimension reduction: PCA
KH (March 2020)
"""
import numpy as np
import pandas as pd
from sklearn.decomposition import KernelPCA

# load gene expression data
dat = pd.read_csv(snakemake.input[0], sep='\t')

# drop gene id column convert to numpy array
dat = dat.iloc[:, :-1].to_numpy()

# perform sample-wise pca projections;
# note that unlike PCA, KernelPCA reduces the column dimension by default
transformer = KernelPCA(n_components=snakemake.config['dimension_reduction']['num_dims'],
                        kernel=snakemake.params['kernel'])
fit = transformer.fit_transform(dat.T)

# convert result to a dataframe with PCs as row and samples as columns
res = pd.DataFrame(fit.T)

res.columns = dat.columns[:-1]
res.index = ['PC' + str(i) for i in range(1, res.shape[0] + 1)]

# store projected data and variance explained
np.savetxt(snakemake.output['fit'], fit, delimiter='\t')
