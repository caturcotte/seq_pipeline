import duckdb


if __name__ == "__main__":
    ref = snakemake.input['ref']
    smp = snakemake.input['smp']
    out_file = snakemake.output[0]

    con = duckdb.connect(snakemake.output['db'], config = {
        'threads': snakemake.threads,
        'memory_limit': (
            f"{snakemake.resources['mem_mb']/1000 - 2}G"
            )
        })
    ref_pqt = con.read_parquet(ref)
    smp_pqt = con.read_parquet(smp)
    con.sql(
    """CREATE
    OR REPLACE VIEW ref AS
    SELECT
        *
    FROM
        ref_pqt
    """
    )
    con.sql(
        """CREATE
        OR REPLACE VIEW temp_vcfs AS
        SELECT
            *
        FROM
            smp_pqt"""
    )
    con.sql(
        """CREATE
        OR REPLACE VIEW vcfs AS
        SELECT
            sample,
            chromosome,
            CAST(position AS INTEGER) AS int_pos,
            reference,
            variant,
            quality,
            genotype,
            depth,
            allele_depth,
            (
                CASE
                    WHEN temp_vcfs.variant = '.' THEN UPPER(temp_vcfs.reference)
                    ELSE UPPER(temp_vcfs.variant)
                END
            ) AS new_variant
        FROM
            temp_vcfs"""
    )
    con.sql(
    """CREATE
    OR REPLACE VIEW samples_rearranged AS
    SELECT
        *,
        (
            CASE
                WHEN ref_allele = new_variant THEN '.'
                ELSE new_variant
            END
        ) AS var_adjusted
    FROM
        vcfs
        INNER JOIN refs ON vcfs.chromosome = refs.chromosome
        AND vcfs.int_pos = refs.int_pos
    WHERE
        new_variant = ref_allele
        OR new_variant = alt_allele"""
    )
    con.sql("SET reserve_insertion_order = true")
    final_pqt = con.sql(
        """CREATE
        OR REPLACE TABLE final_samples AS
        SELECT 
            string_split(sample, '-')[1] AS condition,
            string_split(sample, '-')[2] AS sample_type,
            string_split(sample, '-')[3] AS sample_num,
            UPPER(ref_allele) AS reference,
            UPPER(var_adjusted) AS variant,
            chromosome,
            int_pos AS position,
            CAST(string_split(allele_depth, ',')[1] AS INTEGER) AS ref_reads,
            CAST(string_split(allele_depth, ',')[2] AS INTEGER) AS variant_reads,
            quality AS QUAL,
            (CASE WHEN var_adjusted='.' AND genotype='1/1' THEN '0/0' WHEN var_adjusted !='.' AND genotype='0/0' THEN '1/1' ELSE genotype END) AS GT,
            depth AS DP
        FROM
            samples_rearranged
        ORDER BY
            sample_num, chromosome, position"""
        )
    final_pqt.write_parquet(out_file) 
    con.close()
    # con.sql(
    #     """COPY (
    #         SELECT
    #             *
    #         FROM
    #             final_samples
    #         )
    #     TO
    #         sample_out_file
    #     (FORMAT 'parquet')"""
    # )