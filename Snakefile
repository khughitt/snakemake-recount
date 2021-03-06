"""
Recount2 Pipeline
"""
import os
import pandas as pd
from os.path import join

# get a list of project ids
infile = join(config['output_dir'], 'metadata', 'recount_num_samples.tsv')

if not os.path.exists(infile):
    raise(Exception("Recount Metadata not available! Use the 'fetch_recount_metadata.R' "
                    "script to retrieve the necessary metadata prior to running the pipeline."))

# load project metadata and determine which samples to include
metadata = pd.read_csv(infile, sep = '\t')
project_ids = metadata.project[metadata.number_samples >= config['min_samples']].tolist()

# TESTING
project_ids = project_ids[200:205]

print("Retrieving data for {}/{} projects with sufficient samples.".format(
    len(project_ids), metadata.shape[0]))

rule all:
    input:
        expand(join(config['output_dir'], 'deseq2/{project}_deseq2_results.tsv.gz'), project=project_ids)

rule run_deseq2:
    input:
        counts=join(config['output_dir'], 'datasets/{project}/counts_gene_clean.tsv.gz'),
        clusters=join(config['output_dir'], 'clusters/{project}_clusters.tsv.gz')
    output:
        join(config['output_dir'], 'deseq2/{project}_deseq2_results.tsv.gz')
    run:
        'src/run_deseq2.R'

rule cluster_samples:
    input:
        pca=join(config['output_dir'], 'projected/{project}/pca.tsv.gz'),
        pca_var=join(config['output_dir'], 'projected/{project}/pca_var_explained.tsv'),
        ica=join(config['output_dir'], 'projected/{project}/ica.tsv.gz'),
        kpca_poly=join(config['output_dir'], 'projected/{project}/kpca_poly.tsv.gz'),
        kpca_rbf=join(config['output_dir'], 'projected/{project}/kpca_rbf.tsv.gz')
    output:
        join(config['output_dir'], 'clusters/{project}_clusters.tsv.gz')
    script:
        'src/cluster_samples.py'

rule kpca_rbf:
    input:
        join(config['output_dir'], 'datasets/{project}/counts_gene_clean.tsv.gz')
    output:
        join(config['output_dir'], 'projected/{project}/kpca_rbf.tsv.gz')
    params:
        kernel='rbf'
    script:
        'src/project_kpca.py'

rule kpca_poly:
    input:
        join(config['output_dir'], 'datasets/{project}/counts_gene_clean.tsv.gz')
    output:
        join(config['output_dir'], 'projected/{project}/kpca_poly.tsv.gz')
    params:
        kernel='poly'
    script:
        'src/project_kpca.py'

rule ica:
    input:
        join(config['output_dir'], 'datasets/{project}/counts_gene_clean.tsv.gz')
    output:
        join(config['output_dir'], 'projected/{project}/ica.tsv.gz')
    script:
        'src/project_ica.py'

rule pca:
    input:
        join(config['output_dir'], 'datasets/{project}/counts_gene_clean.tsv.gz')
    output:
        fit=join(config['output_dir'], 'projected/{project}/pca.tsv.gz'),
        var=join(config['output_dir'], 'projected/{project}/pca_var_explained.tsv')
    script:
        'src/project_pca.py'

rule clean_dataset:
    input:
        join(config['output_dir'], 'datasets/{project}/counts_gene.tsv.gz')
    output:
        join(config['output_dir'], 'datasets/{project}/counts_gene_clean.tsv.gz')
    script:
        'src/clean_dataset.R'

rule download_dataset:
    output: 
        counts=join(config['output_dir'], 'datasets/{project}/counts_gene.tsv.gz'),
        phenotype=join(config['output_dir'], 'datasets/{project}/{project}.tsv')
    script:
        'src/download_dataset.R'

