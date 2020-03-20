#!/bin/env Rscript
################################################################################
#
# Performs Differential Expression Analysis using DESeq2
# KH (Mar 2020)
#
################################################################################
library(DESeq2)
library(tidyverse)

# load gene expression data
gene_expr <- read_tsv(snakemake@input$counts, col_types = cols())

gene_mat <- as.matrix(gene_expr[, -ncol(gene_expr)])
rownames(gene_mat) <- gene_expr$gene_id

# load sample cluster data
clusters <- read_tsv(snakemake@input$clusters, col_types = cols())

# method/dim columns
cluster_id_cols <- 1:2

# sanity check
if (!all(colnames(gene_mat) == colnames(clusters)[-cluster_id_cols])) {
  stop("Column mismatch for %s", snakemake@wildcards$project_id)
}

# drop any redundant clusterings
clusters <- clusters[!duplicated(clusters[, -cluster_id_cols]), , drop = FALSE]

# move cluster dimension reduction fields to rownames and transpose so that it is
# in the expected order
clusters <- as.data.frame(clusters)
rownames(clusters) <- sprintf("%s_%s", clusters$method, clusters$dim)
clusters <- clusters[, -cluster_id_cols]

clusters <- as.data.frame(t(clusters))

# convert binary cluster memberships to factors
for (i in 1:ncol(clusters)) {
  clusters[, i] <- factor(clusters[, i], labels = c('A', 'B'))
}

# iterate over clusterings and perform differential expression for each
for (clustering in colnames(clusters)) {
  dds <- DESeqDataSetFromMatrix(gene_mat, clusters, as.formula(sprintf("~%s", clustering)))
  dds <- DESeq(dds)

  # first column correspond to the model intercept; second should be the cluster
  # comparison
  contrast <- resultsNames(dds)[2]

  # sanity check
  if (contrast == 'Intercept') {
    stop('Wrong contrast selected!')
  }

  res <- results(dds, name = contrast) %>%
    as.data.frame() %>%
    rownames_to_column(gene_id) %>%
    select(gene_id, everything())

  write_tsv(res, snakemake@output[[1]])
}

