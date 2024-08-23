import duckdb


if __name__ == "__main__":
    # tuples to avoid SQL injection!
    ref_parent_file = snakemake.input['ref']
    alt_parent_file = snakemake.input['alt']
    ref_parent_name = snakemake.config['ref_parent']
    alt_parent_name = snakemake.config['alt_parent']
    ref_qual_cutoff = snakemake.params['min_qual']
    ref_depth_cutoff = snakemake.params['min_depth']
    ref_out_file = snakemake.output[0]

    con = duckdb.connect(snakemake.output['db'], config = {'threads': snakemake.threads, 'memory_limit': (f"{snakemake.resources['mem_mb']/1000 - 2}G")})
    con.sql("SET preserve_insertion_order = false")
    con.execute("""
        SELECT
            sample,
            chromosome,
            CAST(position AS INTEGER) AS int_pos,
            quality,
            genotype,
            depth,
            allele_depth,
            UPPER(reference) AS up_reference,
            UPPER(variant) AS up_variant,
            (CASE WHEN variant='.' THEN reference ELSE UPPER(variant) END) AS mod_variant
        FROM
            read_parquet([?, ?])""", [ref_parent_file, alt_parent_file])
    parents = con.fetchall()
    con.sql(
        """CREATE
        OR REPLACE VIEW check_unique_variants AS
        SELECT
            chromosome,
            int_pos,
            COUNT(up_variant) AS n_variants
        FROM
            parents
        GROUP BY
            chromosome,
            int_pos
        HAVING
            n_variants >= 2"""
        )
    con.execute("""
        SELECT
            sample,
            parents.chromosome,
            parents.int_pos,
            up_reference AS reference,
            mod_variant,
            quality,
            genotype,
            depth,
            allele_depth
        FROM
            parents
            INNER JOIN check_unique_variants ON (
                parents.chromosome = check_unique_variants.chromosome
                AND parents.int_pos = check_unique_variants.int_pos
            )
        WHERE
            genotype != '0/1'
            AND quality > ?
            AND LENGTH (mod_variant) <= 1
            AND LENGTH (up_reference) <= 1
            AND depth > ?""", [ref_qual_cutoff, ref_depth_cutoff])
    temp_ref = con.fetchall()
    con.execute("""
        SELECT
            sample,
            chromosome,
            int_pos,
            mod_variant AS ref_allele
        FROM
            temp_ref
        WHERE
            sample = ?""", [ref_parent_name])
    ref_parent = con.fetchall()
    con.execute("""
        SELECT
            sample,
            chromosome,
            int_pos,
            mod_variant AS alt_allele
        FROM
            temp_ref
        WHERE
            sample = ?""", [alt_parent_name])
    alt_parent = con.fetchall()
    con.sql(
        """CREATE
        OR REPLACE TABLE ref AS
        SELECT
            *
        FROM
            ref_parent
            INNER JOIN alt_parent ON (
                ref_parent.chromosome = alt_parent.chromosome
                AND ref_parent.int_pos = alt_parent.int_pos
            )"""
        )
    con.sql(
        """CREATE
        OR REPLACE TABLE refs_unique AS
        SELECT
            *
        FROM
            ref
        WHERE
            ref_allele != alt_allele"""
        )
    con.execute(
        """
        SELECT
            chromosome,
            int_pos AS position,
            ref_allele AS reference,
            alt_allele AS variant
        FROM
                ref
        """
    )
    ref_final = con.fetchall()
    ref_final.write_parquet(ref_out_file)
    # con.execute("""COPY (
            # SELECT
                # chromosome,
                # int_pos AS position,
                # ref_allele AS reference,
                # alt_allele AS variant
            # FROM
                # ref
            # )
            # TO
                # ?
            # (FORMAT 'parquet')""", [ref_out_file])
    con.close()