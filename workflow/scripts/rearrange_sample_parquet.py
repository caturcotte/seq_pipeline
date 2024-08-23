#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import argparse
import os

parser = argparse.ArgumentParser(description='Import args from snakemake into iPython')

parser.add_argument('sample_name')
parser.add_argument('ref_file')
parser.add_argument('sample_file')
parser.add_argument('sample_out_file')

args = parser.parse_args()

sample_name = args.sample_name
ref_file = args.ref_file
sample_file = args.sample_file
sample_out_file = args.sample_out_file
os.environ['DATABASE_URL'] = f"duckdb://default.db"

import duckdb
import pandas as pd
import ipywidgets

conn = duckdb.connect(f'{sample_name}.db')

get_ipython().run_line_magic('load_ext', 'sql')
get_ipython().run_line_magic('sql', 'conn --alias duckdb')


# In[ ]:


get_ipython().run_cell_magic('sql', '', '\nSET preserve_insertion_order = false;\n')


# In[ ]:


get_ipython().run_cell_magic('sql', '', "\nCREATE\nOR REPLACE VIEW temp_vcfs AS\nSELECT\n    *\nFROM\n    read_parquet('{{sample_file}}');\n")


# In[ ]:


get_ipython().run_cell_magic('sql', '', "\nCREATE\nOR REPLACE VIEW vcfs AS\nSELECT\n    sample,\n    chromosome,\n    CAST(position AS INTEGER) AS int_pos,\n    reference AS old_ref,\n    variant AS old_variant,\n    quality,\n    genotype,\n    depth,\n    allele_depth,\n    (\n        CASE\n            WHEN temp_vcfs.variant = '.' THEN UPPER(temp_vcfs.reference)\n            ELSE UPPER(temp_vcfs.variant)\n        END\n    ) AS new_variant\nFROM\n    temp_vcfs;\n")


# In[ ]:


get_ipython().run_cell_magic('sql', '', "CREATE\nOR REPLACE TABLE refs AS\nSELECT\n    *\nFROM\n    read_parquet('{{ref_file}}');\n")


# In[ ]:


get_ipython().run_cell_magic('sql', '', "\nCREATE\nOR REPLACE VIEW samples_rearranged AS\nSELECT\n    *,\n    (\n        CASE\n            WHEN reference = new_variant THEN '.'\n            ELSE new_variant\n        END\n    ) AS var_adjusted\nFROM\n    vcfs\n    INNER JOIN refs ON vcfs.chromosome = refs.chromosome\n    AND vcfs.int_pos = refs.position\nWHERE\n    new_variant = reference\n    OR new_variant = variant;\n")


# In[ ]:


get_ipython().run_cell_magic('sql', '', "SET preserve_insertion_order = true;\nCREATE\nOR REPLACE TABLE final_samples AS\nSELECT \n    string_split(sample, '-')[1] AS condition,\n    string_split(sample, '-')[2] AS sample_type,\n    string_split(sample, '-')[3] AS sample_num,\n    UPPER(reference) AS reference,\n    UPPER(var_adjusted) AS variant,\n    chromosome,\n    int_pos AS position,\n    CAST(string_split(allele_depth, ',')[1] AS INTEGER) AS ref_reads,\n    CAST(string_split(allele_depth, ',')[2] AS INTEGER) AS variant_reads,\n    quality AS QUAL,\n    (CASE WHEN var_adjusted='.' AND genotype='1/1' THEN '0/0' WHEN var_adjusted !='.' AND genotype='0/0' THEN '1/1' ELSE genotype END) AS GT,\n    depth AS DP\nFROM\n    samples_rearranged\nORDER BY\n    sample_num, chromosome, position;\n")


# In[ ]:


get_ipython().run_cell_magic('sql', '', "\nCOPY (SELECT * FROM final_samples)\nTO '{{sample_out_file}}'\n(FORMAT 'parquet');\n")

