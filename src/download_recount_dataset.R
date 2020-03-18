#!/bin/env Rscript
#
# Downloads the counts and metadata for a single Recount2 dataset
#
suppressMessages(library(recount))

output_dir <- dirname(snakemake@output$counts)

# download counts and metadata
url <- download_study(snakemake@wildcards$project, type = 'counts-gene', outdir = output_dir)
url <- download_study(snakemake@wildcards$project, type = 'phenotype', outdir = output_dir)

# rename phenotype file
# file.rename(file.path(output_dir, sprintf("%s.tsv", snakemake@wildcards$project)),
#             snakemake@output$phenotype)
