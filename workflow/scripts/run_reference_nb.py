import papermill as pm
import duckdb
import ipywidgets
import os

pm.execute_notebook(
    snakemake.input['old_nb'],
    snakemake.output['new_nb'],
    parameters = dict(
        db=snakemake.output['db'],
        ref_parent_file=snakemake.input['ref'],
        alt_parent_file=snakemake.input['alt'],
        ref_parent_name=snakemake.config['ref_parent'],
        ref_qual_cutoff=snakemake.config['min_qual'],
        ref_depth_cutoff=snakemake.config['min_depth'],
        ref_out_file=snakemake.output[0],
        mem_limit=(snakemake.resources['mem_mb']/1000 - 2)
        )
)