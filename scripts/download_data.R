#!/bin/env Rscript
#
# Downloads the counts and metadata for a single Recount2 dataset
#
suppressMessages(library(recount))

output_dir <- dirname(snakemake@output[[1]])

# download counts and metadata
url <- download_study(snakemake@wildcards$project_id, type = 'counts-gene', outdir = output_dir)
url <- download_study(snakemake@wildcards$project_id, type = 'phenotype', outdir = output_dir)
