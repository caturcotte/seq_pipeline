from pyspark.sql import SparkSession

if __name__ == "__main__":
    spark = SparkSession.builder.appName("VCF DB").getOrCreate()

    ref_parent_file = snakemake.inputs['ref_parent']
    alt_parent_file = snakemake.inputs['alt_parent']
    ref_parent_name = "w1118"
    alt_parent_name = "oregonr"
    ref_qual_cutoff = snakemake.params['qual_cutoff']
    ref_depth_cutoff = snakemake.params['depth_cutoff']
    ref_out_file = snakemake.outputs[0]

    df = spark.read.parquet(ref_parent_file, alt_parent_file)

    df.createOrReplaceTempView("df")

    parents = spark.sql(
        "SELECT sample, chromosome, CAST(position AS INTEGER) AS int_pos, quality, genotype, depth, allele_depth, reference, variant, (CASE WHEN variant='.' THEN reference ELSE variant END) AS mod_variant FROM df"
    )
    parents.createOrReplaceTempView("parents")

    check_unique_variants = spark.sql(
        """SELECT chromosome, int_pos, COUNT(variant) AS n_variants FROM parents
            GROUP BY chromosome, int_pos
            HAVING n_variants >= 2"""
    )
    check_unique_variants.createOrReplaceTempView("check_unique_variants")

    temp_ref_query = """SELECT
        sample,
        parents.chromosome,
        parents.int_pos,
        reference,
        mod_variant,
        quality,
        genotype,
        depth,
        allele_depth
        FROM parents
        INNER JOIN check_unique_variants ON (
            parents.chromosome = check_unique_variants.chromosome
            AND parents.int_pos = check_unique_variants.int_pos
            )
        WHERE genotype!='0/1'
            AND quality > {ref_qual_cutoff}
            AND LENGTH(mod_variant) <= 1
            AND LENGTH(reference) <= 1
            AND depth > {ref_depth_cutoff}"""

    temp_ref = spark.sql(
        temp_ref_query,
        ref_qual_cutoff=ref_qual_cutoff,
        ref_depth_cutoff=ref_depth_cutoff
        )
    temp_ref.createOrReplaceTempView("temp_ref")

    parent_query = """SELECT
    sample,
    chromosome,
    int_pos,
    mod_variant AS {allele}
    FROM temp_ref
    WHERE sample={parent_name}"""

    ref_parent = spark.sql(
        ref_query,
        allele="ref_allele",
        parent_name=ref_parent_name
    )
    ref_parent.createOrReplaceTempView("ref_parent")

    alt_parent = spark.sql(
        ref_query,
        allele="alt_allele",
        parent_name=alt_parent_name
    )
    alt_parent.createOrReplaceTempView("alt_parent")

    ref = spark.sql("""SELECT * FROM ref_parent
        FULL JOIN alt_parent ON (
            ref_parent.chromosome = alt_parent.chromosome
            AND ref_parent.int_pos = alt_parent.int_pos
            )
        ORDER BY ref_parent.chromosome, ref_parent.int_pos"""
    )
    ref.createOrReplaceTempView("ref")

    ref_out = spark.sql(
        """SELECT chromosome, int_pos AS position, ref_allele AS reference, alt_allele AS variant FROM ref""",
        ref=ref
    )
    ref_out.write.parquet(ref_out_file)
