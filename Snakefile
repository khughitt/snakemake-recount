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
        expand(join(config['output_dir'], 'projected/{project}/pca.tsv.gz'), project=project_ids),
        expand(join(config['output_dir'], 'projected/{project}/kpca_linear.tsv.gz'), project=project_ids),
        expand(join(config['output_dir'], 'projected/{project}/kpca_rbf.tsv.gz'), project=project_ids)

rule kpca_rbf:
    input:
        join(config['output_dir'], 'datasets/{project}/counts_gene.tsv.gz')
    output:
        fit=join(config['output_dir'], 'projected/{project}/kpca_rbf.tsv.gz')
    params:
        kernel='rbf'
    script:
        'src/project_kpca.py'

rule kpca_linear:
    input:
        join(config['output_dir'], 'datasets/{project}/counts_gene.tsv.gz')
    output:
        fit=join(config['output_dir'], 'projected/{project}/kpca_linear.tsv.gz')
    params:
        kernel='linear'
    script:
        'src/project_kpca.py'

rule pca:
    input:
        join(config['output_dir'], 'datasets/{project}/counts_gene.tsv.gz')
    output:
        fit=join(config['output_dir'], 'projected/{project}/pca.tsv.gz'),
        var=join(config['output_dir'], 'projected/{project}/pca_var_explained.tsv')
    script:
        'src/project_pca.py'

rule download_recount_project:
    output: 
        counts=join(config['output_dir'], 'datasets/{project}/counts_gene.tsv.gz'),
        phenotype=join(config['output_dir'], 'datasets/{project}/{project}.tsv')
    script:
        'src/download_recount_dataset.R'

