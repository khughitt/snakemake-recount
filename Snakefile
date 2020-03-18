"""
Recount2 Pipeline
"""
import os
import pandas as pd
from os.path import join

# get a list of project ids
infile = join(config['output_dir'], 'metadata', 'recount_url.tsv')

if not os.path.exists(infile):
    raise(Exception("Recount Metadata not available! Use the 'fetch_recount_metadata.R' "
                    "script to retrieve the necessary metadata prior to running the pipeline."))

# load project metadata
metadata = pd.read_csv(infile, sep = '\t')

# get a list of all available project ids
sample_counts = metadata.groupby('project')['path'].count()
project_ids = list(sample_counts.index[sample_counts >= config['min_samples']])

print("Retrieving data for {}/{} projects with sufficient samples.".format(
    len(project_ids), len(sample_counts)))

rule download_all:
    input:
        counts=expand(join(config['output_dir'], 'datasets/{project}/counts_gene.tsv.gz'), project=project_ids),
        phenotype=expand(join(config['output_dir'], 'datasets/{project}/{project}.tsv'), project=project_ids)

rule download_recount_project:
    output: 
        counts=join(config['output_dir'], 'datasets/{project}/counts_gene.tsv.gz'),
        phenotype=join(config['output_dir'], 'datasets/{project}/{project}.tsv')
    script:
        'src/download_recount_dataset.R'

