#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import argparse

parser = argparse.ArgumentParser()
parser.add_argument('sample_file_list')
parser.add_argument('output_file_list')
parser.add_argument('chr_names', nargs="+")


args = parser.parse_args()
sample_file_list = args.sample_file_list
output_file_list = args.output_file_list
chr_names = sorted(args.chr_names)

with open(sample_file_list, 'r') as file:
    sample_files = [i.rstrip() for i in file]
    
with open(output_file_list, 'r') as file:
    output_files = [i.rstrip() for i in file]
print(sample_files)

# In[57]:


import duckdb


# In[ ]:


import pandas as pd


# In[ ]:


import ipywidgets


# In[ ]:


from pathlib import Path


# In[58]:


get_ipython().run_line_magic('reload_ext', 'sql')


# In[59]:


get_ipython().run_line_magic('config', 'SqlMagic.autopandas = False')
get_ipython().run_line_magic('config', 'SqlMagic.feedback = False')
get_ipython().run_line_magic('config', 'SqlMagic.displaycon = False')


# In[60]:


conn = duckdb.connect("reference.db")

get_ipython().run_line_magic('sql', 'conn')


# In[89]:

prelim_inputs = conn.read_parquet(sample_files).df()
prelim_inputs.fillna(0)
prelim_inputs.sort_values(by=['sample', 'chromosome', 'position'])
prelim_inputs = prelim_inputs.loc[prelim_inputs['chromosome'].isin(chr_names)]
prelim_inputs = prelim_inputs.loc[prelim_inputs['DP'] != '.']
prelim_inputs = prelim_inputs.loc[prelim_inputs['QUAL'] != '.']
prelim_inputs = prelim_inputs.loc[prelim_inputs['GT'] != '.']
prelim_inputs['DP'] = prelim_inputs['DP'].astype(int)
prelim_inputs['QUAL'] = prelim_inputs['QUAL'].astype(float)
prelim_inputs['ref_reads'] = prelim_inputs['ref_reads'].astype(int)
prelim_inputs['variant_reads'] = prelim_inputs['variant_reads'].astype(int)

# prelim_inputs.to_csv('test00.csv')
tiger_chrom_name_dict = dict(zip(chr_names, range(1, len(chr_names)+1)))
# import csv
#
# with open("tiger_chrom_dict.csv", "w", newline="") as f:
#     w = csv.DictWriter(f, tiger_chrom_name_dict.keys())
#     w.writeheader()
#     w.writerow(tiger_chrom_name_dict)

# pd.DataFrame(tiger_chrom_name_dict).to_csv('tiger_chrom_dict.csv')
# prelim_inputs = prelim_inputs.astype({
#         'chromosome': 'str',
#         'position': 'int64',
#         'ref_reads': 'int32',
#         'variant_reads': 'int32'
#     })


# Exclude any SNPs where the ref or variant aren't actually in the progeny or are very low coverage over all samples. We want to make sure there's a heterozygous population for every SNP analyzed.

# In[90]:

get_ipython().run_cell_magic('sql', '', 'SET preserve_insertion_order = true;')
get_ipython().run_cell_magic('sql', '', 'CREATE VIEW IF NOT EXISTS good_variants AS\nSELECT\n  chromosome,\n  POSITION,\n  AVG(variant_reads) AS avg_variant_reads,\n  AVG(ref_reads) AS avg_ref_reads\nFROM\n  prelim_inputs\nGROUP BY\n  chromosome,\n  POSITION\nHAVING\n  (\n    avg_variant_reads > 5\n    AND avg_ref_reads > 5\n  ) ORDER BY (chromosome, position);\n')


# In[91]:


# get_ipython().run_cell_magic('sql', '', 'CREATE OR REPLACE TABLE inputs AS\nSELECT * FROM prelim_inputs\nINNER JOIN good_variants ON\n(prelim_inputs.chromosome = good_variants.chromosome AND prelim_inputs.position = good_variants.position);\n')
get_ipython().run_cell_magic('sql', '', 'CREATE OR REPLACE TABLE inputs AS\nSELECT * FROM prelim_inputs\nINNER JOIN good_variants ON\n(prelim_inputs.chromosome = good_variants.chromosome AND prelim_inputs.position = good_variants.position);\n')


get_ipython().run_cell_magic('sql', '', 'inputs << SELECT * FROM inputs;')

# In[105]:

final_inputs = inputs.DataFrame()

final_inputs = final_inputs.sort_values(by=['sample', 'chromosome', 'position'])

# In[109]:


# final_inputs['sample'] = final_inputs['condition'] + '-' + final_inputs['sample_type'] + '-' + final_inputs['sample_num']


# In[110]:


samples = [Path(i).stem.rstrip("_rearranged") for i in sample_files]
tiger_chrom_name_dict = dict(zip(chr_names, range(1, len(chr_names)+1)))
for enum, i in enumerate(samples):
    df_subset = final_inputs.loc[final_inputs['sample'] == i].sort_values(by=['chromosome', 'position'])
    df_subset = df_subset.fillna(0)
    # df_subset = df_subset.loc[df_subset['chromosome'].isin(chr_names)]
    # df_subset = df_subset.replace({'chromosome': tiger_chrom_name_dict})
    # df_subset = df_subset.rename(columns={'chromosome': 'tiger_chrom'})
    df_subset = df_subset.astype({
        'chromosome': 'str',
        'position': 'int64',
        'ref_reads': 'int32',
        'variant_reads': 'int32'
    })
    df_subset = df_subset[['chromosome', 'position', 'reference', 'ref_reads', 'variant', 'variant_reads', 'phase_set', 'QUAL', 'GT']].sort_values(by=['chromosome', 'position'])
    df_subset.to_csv(
        output_files[enum],
        header=True,
        index=False,
        columns=[
            "chromosome",
            "position",
            "reference",
            "ref_reads",
            "variant",
            "variant_reads",
            "phase_set",
            "QUAL",
            "GT"
        ],
    )

