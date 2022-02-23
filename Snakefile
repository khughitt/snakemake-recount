"""
Recount Pipeline
"""
import os
import pandas as pd
from os.path import join

configfile: "config/config.yml"

# get a list of project ids
infile = "metadata/recount_num_samples.tsv"

# load project metadata and determine which samples to include
metadata = pd.read_csv(infile, sep = "\t")
project_ids = metadata.project[metadata.number_samples >= config["min_samples"]].tolist()

# after loading the metadata, set workdir to the desired output directory
workdir: config['out_dir']

projection_methods = ['pca', 'kpca', 'ica']

# TESTING
project_ids = project_ids[200:205]

rule all:
    input: expand('{project_id}/deseq/results.feather', project_id=project_ids)

rule run_deseq2:
    input:
       '{project_id}/filtered/counts_gene.feather',
       '{project_id}/clusters/combined.feather',
    output:
        '{project_id}/deseq/results.feather'
    script:
        'scripts/run_deseq2.R'

rule combine_clustering_results:
    input:
        '{project_id}/clusters/orig.feather',
        expand('{{project_id}}/clusters/{projection}.feather', projection=projection_methods)
    output:
        '{project_id}/clusters/combined.feather'
    script:
        'scripts/combine_clustering_results.py'

rule cluster_projected_samples:
    input:
        '{project_id}/projected/{projection}.feather'
    output:
        '{project_id}/clusters/{projection}.feather'
    script:
        'scripts/cluster_projected_samples.py'

rule cluster_samples:
    input:
       '{project_id}/filtered/counts_gene.feather' 
    output:
       '{project_id}/clusters/orig.feather'
    script:
        'scripts/cluster_samples.py'

rule ica:
    input:
       '{project_id}/filtered/counts_gene.feather' 
    output:
        '{project_id}/projected/ica.feather',
    script:
        'scripts/project_ica.py'

rule kpca:
    input:
       '{project_id}/filtered/counts_gene.feather' 
    output:
        '{project_id}/projected/kpca.feather',
    params:
        kernel='poly'
    script:
        'scripts/project_kpca.py'

rule pca:
    input:
       '{project_id}/filtered/counts_gene.feather' 
    output:
        '{project_id}/projected/pca.feather',
        '{project_id}/projected/pca_var_explained.csv'
    script:
        'scripts/project_pca.py'

rule filter_data:
    input:
       '{project_id}/raw/counts_gene.tsv.gz' 
    output:
       '{project_id}/filtered/counts_gene.feather' 
    script:
        'scripts/filter_data.R'

rule download_data:
    output: 
        '{project_id}/raw/counts_gene.tsv.gz',
        '{project_id}/raw/{project_id}.tsv'
    script:
        'scripts/download_data.R'
