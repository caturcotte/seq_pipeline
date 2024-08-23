import duckdb
import modin.pandas as pd


if __name__ == "__main__":
    tiger_chrom_name_dict = dict(
        zip(
            snakemake.params['chrom_names'],
            range(len(snakemake.params['chrom_names']))
            )
        )

    with duckdb.connect(
        snakemake.output['db'],
        config = {
            'threads': snakemake.threads,
            'memory_limit': (
                f"{snakemake.resources['mem_mb']/1000 - 2}G"
                )
            }
        ) as con:
        df = con.read_parquet(snakemake.input)
        con.sql(
        """COPY 
        """
        )
        df = con.read_parquet(snakemake.input).df()
        df['tiger_chrom'] = df.apply(
            lambda row: tiger_chrom_name_dict[row['chromosome']], axis=1
        )
        df_sql = con.from_df(df)
        con.sql(
        """CREATE OR REPLACE TABLE df_reworked AS
        SELECT
            tiger_chrom,
            position,
            reference,
            ref_reads,
            variant,
            variant_reads
        FROM
            df_sql
        """
        )
        con.sql(
            """COPY 
                (SELECT * FROM df_reworked)
            TO
                '{{snakemake.output[0]}}
            (DELIMITER '\t', HEADER false)"""
        )