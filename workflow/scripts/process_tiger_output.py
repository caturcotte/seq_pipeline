import duckdb
import modin.pandas as pd


def read_dfs(files, db):
    with duckdb.connect(db) as con:
        df = con.read_csv(
            files,
            sep="\t",
            columns = {
                "sample": 'STRING',
                "chromosome": 'STRING',
                "position": 'UBIGINT',
                "base_geno": 'STRING',
                "hmm_state1": 'STRING',
                "hmm_state2": 'STRING',
                "reference": 'STRING',
                "ref_reads": 'UINTEGER',
                "variant": 'STRING',
                "var_reads": 'UINTEGER'
            }
        ).df()
    return df


def create_master_df(df, genotype_dict):
    df = df.replace(genotype_dict)
    #TODO: figure out how to do this for F1 progeny
    #for our specific cross scheme to ID crossovers, heterozygous and homozygous CB4856 calls are should all be considered CB4856
    df['hmm_state1'] = df['hmm_state1'].replace(
        {'het':alt_parental_genotype}
        )
    
    df['chrom_id'] = df.apply(
        lambda row: row['sample']+"-"+str(row['chromosome']),
        axis=1
        )
    return df


def make_transition_intervals_df(master_df):
    
    #These new transition intervals are the potential crossover intervals that can be classified later as crossover or 'not crossover'

    new_dfs = []

    for chrom_id in master_df.chrom_id.unique():
        df = master_df[master_df['chrom_id']==chrom_id]
        starts = list(df.start)
        stops = list(df.stop)
        states = list(df.hmm_state)
        new_intervals = []

        for i in range(len(stops)-1):
            if states[i] != states[i+1]:
                row = [
                    df['sample'].unique()[0],
                    df['chromosome'].unique()[0],
                    stops[i],
                    starts[i+1],
                    'transition',
                    df['chrom_id'].unique()[0],
                    (starts[i+1]-stops[i])
                    ]
                new_intervals.append(row)
        new_intervals_df = pd.DataFrame(new_intervals, columns=df.columns)
        full_intervals_df = pd.concat(
            [df, new_intervals_df]
            ).sort_values(by='start').reset_index(drop=True)
        new_dfs.append(full_intervals_df)
    return pd.concat(new_dfs)


def create_state_intervals_df(df, genotype_dict):
    #This will create a much smaller dataframe than the marker dataframe above by collapsing runs of the same HMM state into intervals and creating a new row for a 2-marker transition interval between HMM states

    #for our specific cross scheme to ID crossovers, heterozygous and homozygous CB4856 calls are should all be considered CB4856


    df['hmm_state'] = df.apply(
        lambda row: genotype_dict[(row['hmm_state'])], axis=1
        )
    df['chrom_id'] = df.apply(
        lambda row: row['sample']+"-"+str(row['chromosome']), axis=1
        )
    df['length'] = df.apply(
        lambda row: row['stop']-row['start']+1, axis=1
        )
    df = make_transition_intervals_df(df)
    df = df.sort_values(
        by=['chrom_id', 'chromosome', 'start']
        ).reset_index(drop=True)
    return df


def get_marker_counts_per_interval(master_bases_df, master_intervals_df):
    
    #This is useful to retain marker density information for each interval, which is needed for later maths and classification of crossover intervals

    new_dfs = []

    for chrom_id in master_bases_df.chrom_id.unique():
        bases_df = master_bases_df[master_bases_df['chrom_id']==chrom_id]
        intervals_df = master_intervals_df[
            master_intervals_df['chrom_id']==chrom_id
            ]
        starts = list(intervals_df.start)
        stops = list(intervals_df.stop)
        counts = []

        for i in range(len(starts)):
            counts.append(len(bases_df[(bases_df['position']>=starts[i]) & (bases_df['position']<=stops[i]) ]))
        intervals_df['marker_counts'] = counts
        new_dfs.append(intervals_df)
    markers_counted_df = pd.concat(new_dfs)
    return markers_counted_df


if __name__ == "__main__":
    ref_parental_genotype = snakemake.config['ref_parent']
    alt_parental_genotype = snakemake.config['alt_parent']
    genotype_dict = {
        'CC':ref_parental_genotype,
        "CL":"het",
        "LL":alt_parental_genotype,
        "CU":"u"+ref_parental_genotype,
        "LU":'u'+alt_parental_genotype,
        "UN":"unknown",
        '?':"?"
        }
    df = read_dfs(snakemake.input, snakemake.output['db'])
    tiger_marker_df = create_TIGER_master_df(
        df, genotype_dict, snakemake.output['db']
        )
    tiger_pre_intervals = create_state_intervals_df(
        tiger_marker_df, genotype_dict
        )
    tiger_intervals = get_marker_counts_per_interval(
        tiger_marker_df, tiger_pre_intervals
        )
    tiger_intervals.to_parquet(snakemake.output['marker_df'])
    tiger_intervals.to_parquet(snakemake.output['pre_intervals_df'])
    tiger_intervals.to_parquet(snakemake.output['intervals_df'])