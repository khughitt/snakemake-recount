#!/bin/env Rscript
#
# Downloads the latest Recount2 metadata
#
suppressMessages(library(recount))
suppressMessages(library(readr))
suppressMessages(library(yaml))

# create output directory
config <- read_yaml('../config/config-v1.0.yml')
output_dir <- file.path(config$output_dir, 'metadata')

if (!dir.exists(output_dir)) {
  dir.create(output_dir, mode = '0755', recursive = TRUE)
}

# retrieve detailed metadata for sra, tcga, gtex
write_tsv(as.data.frame(all_metadata('sra')), file.path(output_dir, 'sra.tsv'))
write_tsv(as.data.frame(all_metadata('gtex')), file.path(output_dir, 'gtex.tsv'))
write_tsv(as.data.frame(all_metadata('tcga')), file.path(output_dir, 'tcga.tsv'))

# retrieve summary metadata including all supported projects
write_tsv(recount_url, file.path(output_dir, 'recount_url.tsv'))

