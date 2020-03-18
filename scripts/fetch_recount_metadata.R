#!/bin/env Rscript
#
# Downloads the latest Recount2 metadata
#
suppressMessages(library(recount))
suppressMessages(library(tidyverse))
suppressMessages(library(yaml))

# create output directory
config <- read_yaml('../config/config-v1.0.yml')
output_dir <- file.path(config$output_dir, 'metadata')

if (!dir.exists(output_dir)) {
  dir.create(output_dir, mode = '0755', recursive = TRUE)
}

# the "recount_abstract" object contain information about the number of
# samples in each study, as well as abstracts, etc.
sra <- as.data.frame(all_metadata('sra'))
gtex <- as.data.frame(all_metadata('gtex'))
tcga <- as.data.frame(all_metadata('tcga'))

sra <- sra %>%
  inner_join(recount_abstract, by = 'project')
gtex <- gtex %>%
  inner_join(recount_abstract, by = 'project')
tcga <- tcga %>%
  inner_join(recount_abstract, by = 'project')

# create a mapping from project ids to number of samples
num_samples <- sra %>%
  select(project, number_samples) %>%
  unique()

num_samples <- rbind(num_samples, c('TCGA', tcga$number_samples[1]))
num_samples <- rbind(num_samples, c('GTEX', gtex$number_samples[1]))

# store detailed metadata for sra, tcga, gtex
write_tsv(sra, file.path(output_dir, 'sra.tsv'))
write_tsv(gtex, file.path(output_dir, 'gtex.tsv'))
write_tsv(tcga, file.path(output_dir, 'tcga.tsv'))

# store dataset-level for all recount2 projects
write_tsv(recount_url, file.path(output_dir, 'recount_url.tsv'))

# store project sample numbers
write_tsv(num_samples, file.path(output_dir, 'recount_num_samples.tsv'))

# store gene annotations
write_tsv(recount_genes, file.path(output_dir, 'recount_genes.tsv'))

