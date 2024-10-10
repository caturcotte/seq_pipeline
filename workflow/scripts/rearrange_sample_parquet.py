#!/usr/bin/env python
# coding: utf-8

# In[1]:


import argparse
import os

parser = argparse.ArgumentParser(description='Import args from snakemake into iPython')

parser.add_argument('sample_name')
parser.add_argument('ref_file')
parser.add_argument('paternal_file')
parser.add_argument('sample_file')
parser.add_argument('sample_out_file')

args = parser.parse_args()

sample_name = args.sample_name
ref_file = args.ref_file
paternal_file = args.paternal_file
sample_file = args.sample_file
sample_out_file = args.sample_out_file

# sample_name = 'WT-F1-018'
# ref_file = '/work/users/c/a/cannecar/co_analysis/seq_pipeline/data/parquets/reference.parquet'
# sample_file = '/work/users/c/a/cannecar/co_analysis/seq_pipeline/data/parquets/WT-F1-018_converted.parquet'
# sample_out_file = 'test_sample.pq'
# os.environ['DATABASE_URL'] = f"duckdb://default.db"


# Sometimes ipython throws a hissy about importing duckdb or loading the database. I broke it up into multiple cells because for some reason that seems to help?

# In[ ]:


import duckdb


# In[ ]:


import pandas as pd


# In[ ]:


import ipywidgets


# In[ ]:


conn = duckdb.connect(f'{sample_name}.db')


# In[ ]:


get_ipython().run_line_magic('reload_ext', 'sql')


# In[ ]:


get_ipython().run_line_magic('sql', 'conn --alias duckdb')


# In[2]:


get_ipython().run_cell_magic('sql', '', 'SET\n  preserve_insertion_order = FALSE;\n')


# Load in the SNPs from the sample file.

# In[3]:


get_ipython().run_cell_magic('sql', '', "CREATE TABLE IF NOT EXISTS temp_vcfs AS\nSELECT\n  *\nFROM\n  read_parquet ('{{sample_file}}');\n")


# Make a new table where '.' in the variants column is replaced with the reference allele. Also rename some of the columns to make it easier to keep track of them later on.

# In[5]:


get_ipython().run_cell_magic('sql', '', "CREATE TABLE IF NOT EXISTS vcfs AS\nSELECT\n  sample,\n  chromosome,\n  CAST(position AS INTEGER) AS int_pos,\n  reference AS old_ref,\n  variant AS old_variant,\n  quality,\n  genotype AS old_genotype,\n  depth,\n  allele_depth,\n  phase_set, (\n    CASE\n      WHEN temp_vcfs.variant = '.' THEN UPPER(temp_vcfs.reference)\n      ELSE UPPER(temp_vcfs.variant)\n    END\n  ) AS new_variant\nFROM\n  temp_vcfs;\n")


# Split the allele depth (format "ref_reads,variant_reads") into ref reads and variant reads.

# In[7]:


get_ipython().run_cell_magic('sql', '', "CREATE TABLE IF NOT EXISTS ad_split AS \nSELECT\n    *,\n    string_split(allele_depth, ',')[1] AS initial_ref_reads,\n    string_split(allele_depth, ',')[2] AS initial_variant_reads,\nFROM vcfs;\n")


# Load in the reference parent parquet file.

# In[12]:


get_ipython().run_cell_magic('sql', '', "CREATE TABLE IF NOT EXISTS refs AS\nSELECT\n  *\nFROM\n  read_parquet ('{{ref_file}}');\n")


# Join the parent and sample files and filter for only variants that either match the reference or variant of the parent.

# In[14]:


get_ipython().run_cell_magic('sql', '', 'CREATE TABLE IF NOT EXISTS samples_rearranged AS\nSELECT\n  *\nFROM\n  ad_split\n  INNER JOIN refs ON ad_split.chromosome = refs.chromosome\n  AND ad_split.int_pos = refs.position\nWHERE\n  new_variant = reference\n  OR new_variant = variant;\n')


# Load the paternal SNPs file into the database and filter out any variants on chromosomes other than chr2 that are found in the paternal genome.

# In[21]:


get_ipython().run_cell_magic('sql', '', 'CREATE TABLE IF NOT EXISTS paternal_snps AS\nSELECT\n    *,\n    UPPER(reference) AS pat_ref,\n    (\n        CASE\n            WHEN variant = \'.\' THEN UPPER(pat_ref)\n            ELSE UPPER(variant) END\n    ) AS pat_var\nFROM\n    read_parquet(\'{{paternal_file}}\');\n')


get_ipython().run_cell_magic('sql', '', "CREATE\nOR REPLACE TABLE paternal_snps_filtered AS\nSELECT\n    *\nFROM samples_rearranged\nLEFT JOIN\n    paternal_snps\nON (\n    samples_rearranged.chromosome = paternal_snps.chromosome\n    AND samples_rearranged.position = paternal_snps.position\n    )\nWHERE\n    samples_rearranged.chromosome LIKE 'chr2%' OR pat_var != new_variant;\n")


# Swap the ref read and variant read counts based on whether the reference parent SNP matched the original reference or not.

# In[20]:


get_ipython().run_cell_magic('sql', '', "CREATE\nOR REPLACE TABLE ads_swapped AS\nSELECT\n  sample,\n  chromosome,\n  position,\n  reference,\n  variant,\n  quality,\n phase_set,  (\n    CASE\n      WHEN old_ref = reference THEN old_genotype\n      WHEN old_ref != reference\n      AND old_genotype = '1/1' THEN '0/0'\n      WHEN old_ref != reference\n      AND old_genotype = '0/0' THEN '1/1'\n      WHEN old_genotype = '0/1' THEN '0/1'\n     WHEN old_genotype = '0|1' AND old_ref != reference THEN '1|0' WHEN old_genotype = '1|0' AND old_ref != reference THEN '0|1' \n      ELSE './.'\n    END\n  ) AS genotype,\n  depth,\n  (\n    CASE\n      WHEN reference = old_ref THEN initial_ref_reads\n      WHEN reference = new_variant THEN initial_variant_reads\n      ELSE '0'\n    END\n  ) AS ref_reads,\n  (\n    CASE\n      WHEN old_ref = reference THEN initial_variant_reads\n      WHEN reference = new_variant THEN initial_ref_reads\n      ELSE '0'\n    END\n  ) AS variant_reads\nFROM\n  paternal_snps_filtered;\n")


# In[22]:


get_ipython().run_cell_magic('sql', '', 'SET preserve_insertion_order = true;\nCREATE\nOR REPLACE TABLE final_samples AS\nSELECT\n    sample,\n    UPPER(reference) AS reference,\n    UPPER(variant) AS variant,\n    chromosome,\n    position,\n    ref_reads,\n    variant_reads,\n    quality AS QUAL,\n    genotype AS GT,\n    depth AS DP, phase_set\nFROM\n    ads_swapped\nORDER BY\n    sample, chromosome, position;\n')


# In[24]:


get_ipython().run_cell_magic('sql', '', "COPY (\n  SELECT\n    *\n  FROM\n    final_samples\n) TO '{{sample_out_file}}' (FORMAT 'parquet');\n")

