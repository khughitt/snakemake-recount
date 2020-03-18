"""
Dimension reduction: ICA
KH (March 2020)
"""
import pandas as pd
from sklearn.decomposition import FastICA

# load gene expression data
dat = pd.read_csv(snakemake.input[0], sep='\t')

sample_ids = dat.columns[:-1]

# drop gene id column convert to numpy array
dat = dat.iloc[:, :-1].to_numpy()

# perform sample-wise ICA projections;

# note that by default, PCA reduces column dimensionality, while KPCA/FastICA reduce row
# dimensionality
transformer = FastICA(n_components=snakemake.config['dimension_reduction']['num_dims'],
                      random_state=0)
fit = transformer.fit_transform(dat.T)

# convert result to a dataframe with PCs as row and samples as columns
res = pd.DataFrame(fit.T)

res.columns = sample_ids
res.index = ['PC' + str(i) for i in range(1, res.shape[0] + 1)]

# store projected data and variance explained
res.to_csv(snakemake.output['fit'], sep='\t')
