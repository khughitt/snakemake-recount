"""
Dimension reduction: PCA
KH (March 2020)
"""
import pandas as pd
from sklearn.cluster import KMeans

# load projected data
pca = pd.read_csv(snakemake.input['pca'], sep='\t', index_col=0)
ica = pd.read_csv(snakemake.input['ica'], sep='\t', index_col=0)
kpca_linear = pd.read_csv(snakemake.input['kpca_linear'], sep='\t', index_col=0)
kpca_rbf = pd.read_csv(snakemake.input['kpca_rbf'], sep='\t', index_col=0)

# initialize KMeans engine
kmeans = KMeans(n_clusters=2, random_state=0)

# result list to store cluster labels
res = []

# iterate over PCs/components and perform 1d kmeans clustering on each
for i in range(pca.shape[0]):
    # PCA
    dat = pca.iloc[i, :].to_numpy().reshape(-1,1)
    res.append(['PCA', pca.index[i]] + list(kmeans.fit(dat).labels_))

    # ICA
    dat = ica.iloc[i, :].to_numpy().reshape(-1,1)
    res.append(['ICA', ica.index[i]] + list(kmeans.fit(dat).labels_))

    # KPCA (linear)
    dat = kpca_linear.iloc[i, :].to_numpy().reshape(-1,1)
    res.append(['KPCA-linear', kpca_linear.index[i]] + list(kmeans.fit(dat).labels_))

    # KPCA (RBF)
    dat = kpca_rbf.iloc[i, :].to_numpy().reshape(-1,1)
    res.append(['KPCA-rbf', kpca_rbf.index[i]] + list(kmeans.fit(dat).labels_))

# create a dataframe of cluster assignments
res = pd.DataFrame(res)
res.columns = ['method', 'dim'] + list(pca.columns)

# store result
res.to_csv(snakemake.output[0], sep='\t', index=False)
