#!/bin/env Rscript
################################################################################
#
# Cleans a Recount2 dataset by removing any outlier samples with very low sample
# correlation or too few reads, and removing genes with zero variance.
#
# KH (Mar 2020)
#
################################################################################
suppressMessages(library(tidyverse))

# load recount2 count data
gene_expr <- read_tsv(snakemake@input[[1]], col_types = cols())

orig_size <- ncol(gene_expr)

# first, remove any genes with zero variance as these are uninformative and may cause
# problems for methods like KPCA
row_vars <- apply(gene_expr[, -ncol(gene_expr)], 1, var)
mask <- row_vars > 0

gene_expr <- gene_expr[mask, ]

# check for samples with too few reads and, if found, remove them;
# here, we define that to be any samples with less than 10% of the median number of
# reads for a sample
sample_counts <- colSums(gene_expr[, -ncol(gene_expr)])

count_ratios <- sample_counts / median(sample_counts) 
mask <- c(count_ratios >= snakemake@config$sample_filtering$min_sample_size_ratio, TRUE)

gene_expr <- gene_expr[, mask]

# next, compute the pairwise sample correlations, and remove
cor_mat <- cor(gene_expr[, -ncol(gene_expr)] * 1.0)

median_pairwise_cor <- apply(cor_mat, 2, median)

# remove samples with average correlations less than <N> standard deviations from the
# mean of the median sample pairwise correlations
cutoff <- mean(median_pairwise_cor) - snakemake@config$sample_filtering$sample_cor_iqr_cutoff * sd(median_pairwise_cor)

# remove outlier samples with low intra-sample correlations
mask <- c(median_pairwise_cor >= cutoff, TRUE)
gene_expr <- gene_expr[, mask]

write_tsv(gene_expr, snakemake@output[[1]])

save.image('~/tmp.rda')

# save filtered version of count data
if (ncol(gene_expr) < orig_size) {
  print(sprintf("[INFO] Removed %d/%d outlier samples for dataset %s", orig_size - ncol(gene_expr),
        orig_size - 1, snakemake@wildcards$project))
}

