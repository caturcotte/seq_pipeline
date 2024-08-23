#!/usr/bin/env python
# coding: utf-8

# Connect to the SQL database.

# In[ ]:


import argparse
import os

parser = argparse.ArgumentParser(description='Import args from snakemake into iPython')

parser.add_argument('ref_parent_file')
parser.add_argument('alt_parent_file')
parser.add_argument('ref_parent_name')
parser.add_argument('alt_parent_name')
parser.add_argument('ref_qual_cutoff')
parser.add_argument('ref_depth_cutoff')
parser.add_argument('ref_out_file')
parser.add_argument('db')

args = parser.parse_args()

ref_parent_file = args.ref_parent_file
alt_parent_file = args.alt_parent_file
ref_qual_cutoff = args.ref_qual_cutoff
ref_depth_cutoff = args.ref_depth_cutoff
ref_out_file = args.ref_out_file
os.environ['DATABASE_URL'] = f"duckdb://default.db"


# In[ ]:


import duckdb
import pandas as pd
import ipywidgets

# %config SqlMagic.autopandas = True
# %config SqlMagic.feedback = False
# %config SqlMagic.displaycon = False

get_ipython().run_line_magic('load_ext', 'sql')

# get_ipython().run_line_magic('config', 'SqlMagic.autopandas = False')
# get_ipython().run_line_magic('config', 'SqlMagic.feedback = False')
# get_ipython().run_line_magic('config', 'SqlMagic.displaycon = False')

conn = duckdb.connect()
get_ipython().run_line_magic('sql', 'conn --alias duckdb')


# In[43]:


get_ipython().run_cell_magic('sql', '', 'SET preserve_insertion_order = false;\n')


# Define the locations of the references and samples, and quality/depth cutoffs.

# In[44]:


# ref_parent_file = snakemake.input['ref']
# alt_parent_file = snakemake.input['alt']
# ref_parent_name = snakemake.config['ref_parent']
# alt_parent_name = snakemake.config['alt_parent']

# ref_qual_cutoff = snakemake.params['min_qual']
# ref_depth_cutoff = snakemake.params['min_depth']

# ref_out_file = snakemake.output[0]


# Import the reference files as a table. 

# In[45]:


get_ipython().run_cell_magic('sql', '', "CREATE\nOR REPLACE VIEW parents AS\nSELECT\n    sample,\n    chromosome,\n    CAST(position AS INTEGER) AS int_pos,\n    quality,\n    genotype,\n    depth,\n    allele_depth,\n    UPPER(reference) AS up_reference,\n    UPPER(variant) AS up_variant,\n    (CASE WHEN variant='.' THEN reference ELSE UPPER(variant) END) AS mod_variant\nFROM\n    read_parquet(['{{ref_parent_file}}', '{{alt_parent_file}}']);\n")


# Isolate sites where there are 2 variants (meaning one parent is different from the other).

# In[46]:


get_ipython().run_cell_magic('sql', '', 'CREATE\nOR REPLACE VIEW check_unique_variants AS\nSELECT\n    chromosome,\n    int_pos,\n    COUNT(up_variant) AS n_variants\nFROM\n    parents\nGROUP BY\n    chromosome,\n    int_pos\nHAVING\n    n_variants >= 2;\n')


# Create new reference from the reference parent.

# In[47]:


get_ipython().run_cell_magic('sql', '', "CREATE\nOR REPLACE VIEW temp_ref AS\nSELECT\n    sample,\n    parents.chromosome,\n    parents.int_pos,\n    up_reference AS reference,\n    mod_variant,\n    quality,\n    genotype,\n    depth,\n    allele_depth\nFROM\n    parents\n    INNER JOIN check_unique_variants ON (\n        parents.chromosome = check_unique_variants.chromosome\n        AND parents.int_pos = check_unique_variants.int_pos\n    )\nWHERE\n    genotype != '0/1'\n    AND quality > '{{ref_qual_cutoff}}'\n    AND LENGTH (mod_variant) <= 1\n    AND LENGTH (up_reference) <= 1\n    AND depth > '{{ref_depth_cutoff}}';\n")


# In[48]:


get_ipython().run_cell_magic('sql', '', "CREATE\nOR REPLACE VIEW ref_parent AS\nSELECT\n    sample,\n    chromosome,\n    int_pos,\n    mod_variant AS ref_allele\nFROM\n    temp_ref\nWHERE\n    sample = '{{ref_parent_name}}';\n\nCREATE\nOR REPLACE VIEW alt_parent AS\nSELECT\n    sample,\n    chromosome,\n    int_pos,\n    mod_variant AS alt_allele\nFROM\n    temp_ref\nWHERE\n    sample = '{{alt_parent_name}}';\n")


# In[49]:


get_ipython().run_cell_magic('sql', '', 'CREATE\nOR REPLACE TABLE ref AS\nSELECT\n    *\nFROM\n    ref_parent\n    INNER JOIN alt_parent ON (\n        ref_parent.chromosome = alt_parent.chromosome\n        AND ref_parent.int_pos = alt_parent.int_pos\n    );\n')


# In[50]:


get_ipython().run_cell_magic('sql', '', 'CREATE\nOR REPLACE TABLE refs_unique AS\nSELECT\n    *\nFROM\n    ref\nWHERE\n    ref_allele != alt_allele;\n')


# In[58]:


get_ipython().run_cell_magic('sql', '', "COPY (SELECT chromosome, int_pos AS position, ref_allele AS reference, alt_allele AS variant FROM refs_unique)\nTO '{{ref_out_file}}'\n(FORMAT 'parquet');\n")

