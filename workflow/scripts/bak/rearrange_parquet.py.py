
import duckdb

ref_parent_file = snakemake.input['ref']
alt_parent_file = snakemake.input['alt']
sample_file_glob = snakemake.input['prog']
ref_parent_name = snakemake.config['ref_parent']
alt_parent_name = snakemake.config['alt_parent']
ref_qual_cutoff = 200
ref_depth_cutoff = 30
ref_out_file = snakemake.output['ref']
sample_out_file = snakemake.output['prog']

with duckdb.connect(
    "vcf_dfs.db", config={'threads': snakemake.threads}
    ) as con:
    con.sql("SET memory_limit = '10GB'")
    con.sql("SET enable_progress_bar = true")
    con.sql("SELECT * FROM parents GROUP BY ANY_VALUE(sample), chromosome, int_pos"),