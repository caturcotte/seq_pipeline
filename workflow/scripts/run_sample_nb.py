import papermill as pm

pm.execute_notebook(
   snakemake.input['old_nb'],
   snakemake.output['new_nb'],
   parameters = dict(
       db=f'duckdb:///{snakemake.output['db']}',
        sample_file = snakemake.input['smp'],
        ref_file = snakemake.input['ref'],
        sample_out_file = snakemake.output[0]
   )
)