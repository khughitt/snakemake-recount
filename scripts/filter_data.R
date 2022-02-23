#!/bin/env Rscript
################################################################################
#
# Cleans a Recount2 dataset by removing any outlier samples with very low sample
# correlation or too few reads, and removing genes with zero variance.
#
# kh (mar.2020)
#
################################################################################
suppressMessages(library(arrow))
suppressMessages(library(tidyverse))

# load recount2 count data
gene_expr <- read_tsv(snakemake@input[[1]], col_types = cols()) %>%
  column_to_rownames('gene_id')

orig_size <- ncol(gene_expr)

# first, remove any genes with zero variance as these are uninformative and may cause
# problems for methods like KPCA
row_vars <- apply(gene_expr, 1, var)
mask <- row_vars > 0

if (sum(!mask) > 0) {
  print(sprintf("Excluding %d/%d genes with 0 variance before sample filtering (%s)",
                sum(!mask), length(mask), snakemake@wildcards$project_id))
  gene_expr <- gene_expr[mask, ]
}

# check for samples with too few reads and, if found, remove them;
# here, we define that to be any samples with less than 10% of the median number of
# reads for a sample
sample_counts <- colSums(gene_expr)

count_ratios <- sample_counts / median(sample_counts)
mask <- count_ratios >= snakemake@config$sample_filtering$min_sample_size_ratio

if (sum(!mask) > 0) {
  print(sprintf("Excluding %d/%d samples with too few reads (%s)",
                sum(!mask), length(mask), snakemake@wildcards$project_id))
  gene_expr <- gene_expr[, mask]
}

# next, compute the pairwise sample correlations, and remove
if (snakemake@config$sample_filtering[['enable_sample_cor_stdev_filter']]) {
  cor_mat <- cor(gene_expr * 1.0)

  median_pairwise_cor <- apply(cor_mat, 2, median)

  # remove samples with average correlations greater than <N> standard deviations from the
  # mean pairwise correlation score across all pairs of samples
  num_sds <- snakemake@config$sample_filtering$sample_cor_stdev_cutoff
  cutoff <- mean(median_pairwise_cor) - num_sds * sd(median_pairwise_cor)

  # remove outlier samples with low intra-sample correlations
  mask <- median_pairwise_cor >= cutoff

  if (sum(!mask) > 0) {
    print(sprintf("Excluding %d/%d outlier samples (%s)",
                  orig_size - ncol(gene_expr), orig_size, snakemake@wildcards$project_id))
    gene_expr <- gene_expr[, mask]
  }
}

# remove genes with too many 0's
num_zeros <- apply(gene_expr, 1, function(x) {
  sum(x == 0)
})
cutoff <- round(snakemake@config$gene_filtering$max_zeros * ncol(gene_expr))

mask <- num_zeros <= cutoff

if (sum(!mask) > 0) {
  print(sprintf("Excluding %d/%d genes with more than %d 0's (%s)",
                sum(!mask), length(mask), cutoff, snakemake@wildcards$project_id))
  gene_expr <- gene_expr[mask, ]
}

# check one more time for any genes that have 0 variance as a result of the sample
# filtering
row_vars <- apply(gene_expr, 1, var)
mask <- row_vars > 0

if (sum(!mask) > 0) {
  print(sprintf("Excluding %d/%d genes with 0 variance after sample filtering (%s).",
                sum(!mask), length(mask), snakemake@wildcards$project_id))
  gene_expr <- gene_expr[mask, ]
}

# save filtered version of count data
gene_expr %>%
  rownames_to_column("ensgene") %>%
  write_feather(snakemake@output[[1]])
