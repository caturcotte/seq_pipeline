import modin.pandas as pd


tiger_chrom_name_dict = dict(
    zip(
        snakemake.params['chrom_names'],
        range(len(snakemake.params['chrom_names']))
        )
    )

dfs = []
for i in snakemake.input:
    df = pd.read_parquet(i)
    df['tiger_chrom'] = df.apply(
        lambda row: tiger_chrom_name_dict[row['chromosome']], axis=1
    )
    dfs.append(df)
all_dfs = pd.concat(dfs)
dfs_final = all_dfs[['tiger_chrom', 'position', 'reference', 'ref_reads', 'variant', 'variant_reads']]
dfs_final.write_csv(snakemake.output[0], sep='\t', header=False)