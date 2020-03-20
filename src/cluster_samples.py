"""
Dimension reduction: PCA
KH (March 2020)
"""
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

# load projected data
proj_dat = {
    "PCA": pd.read_csv(snakemake.input['pca'], sep='\t', index_col=0),
    "ICA": pd.read_csv(snakemake.input['ica'], sep='\t', index_col=0),
    "KPCA-poly": pd.read_csv(snakemake.input['kpca_poly'], sep='\t', index_col=0),
    "KPCA-rbf": pd.read_csv(snakemake.input['kpca_rbf'], sep='\t', index_col=0)
}

# initialize KMeans engine
kmeans = KMeans(n_clusters=2, random_state=0)

# result list to store cluster labels
res = []

# iterate over data projections and PCs/components and perform 1d kmeans clustering on each
num_dims = proj_dat['PCA'].shape[0]

for method in proj_dat:
    for i in range(num_dims):
        dat = proj_dat[method].iloc[i, :].to_numpy().reshape(-1, 1)
        res.append([method, i + 1] + list(kmeans.fit(dat).labels_))

        # TESTING
        #  if len(set(list(kmeans.fit(dat).labels_))) == 1:
        #      sys.exit("Only one cluster found for {} {}".format(method, i))

# create a dataframe of cluster assignments
res = pd.DataFrame(res)
res.columns = ['method', 'dim'] + list(proj_dat['PCA'].columns)

# for each projection, genarate a scatterplot of dim 1 vs. dim 2
for method in proj_dat:
    # get projected data values for first two dimensions
    dim1_values = proj_dat[method].iloc[0, :].tolist()
    dim2_values = proj_dat[method].iloc[1, :].tolist()

    # get sample cluster assignments for the first two projected dimensions
    dim1_clusters = res[(res.method == method) & (res.dim == 1)].iloc[0, 2:].tolist()
    dim2_clusters = res[(res.method == method) & (res.dim == 2)].iloc[0, 2:].tolist()

    plt_dat = pd.DataFrame({
        "dim1": dim1_values,
        "dim2": dim2_values,
        "clusters1": dim1_clusters,
        "clusters2": dim2_clusters
    })

    # determine axis labels to use
    if method == 'ICA':
        xlab = 'Dim 1'
        ylab = 'Dim 2'
    else:
        xlab = 'PC1'
        ylab = 'PC2'

    # get index of method in snakemake wildcards (e.g. "KPCA-poly" -> "kpca_poly")
    ind = method.lower().replace('-', '_')

    # plot output filepath template
    outfile_prefix = snakemake.input[ind].replace('.tsv.gz', '')

    # for regular PCA, we can also include the % variance explained by each dimension
    if method == 'PCA':
        pca_var = pd.read_csv(snakemake.input['pca_var'], sep='\t', index_col=0)

        xlab = "{} (var explained: {:0.3} %)".format(xlab, pca_var.iloc[0, 0] * 100)
        ylab = "{} (var explained: {:0.3} %)".format(ylab, pca_var.iloc[1, 0] * 100)

    # generate scatterplots of first two dimensions with points colored by clusterings
    # based on each of the dimensions

    # Dim1 plot
    ax = sns.scatterplot(x="dim1", y="dim2", data=plt_dat, hue="clusters1")

    plt_title = "{}: {} (Dimension 1 clusters)".format(snakemake.wildcards['project'], method)
    ax.set(xlabel=xlab, ylabel=ylab, title=plt_title)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:3], labels[:3])

    fig = ax.get_figure()
    fig.savefig(outfile_prefix + '_1.png')

    # Dim2 plot
    ax = sns.scatterplot(x="dim1", y="dim2", data=plt_dat, hue="clusters2")

    plt_title = "{}: {} (Dimension 2 clusters)".format(snakemake.wildcards['project'], method)
    ax.set(xlabel=xlab, ylabel=ylab, title=plt_title)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[:3], labels[:3])

    fig = ax.get_figure()
    fig.savefig(outfile_prefix + '_2.png')

# store result
res.to_csv(snakemake.output[0], sep='\t', index=False)
