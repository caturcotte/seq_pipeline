SELECT getenv('REF_PARENT_FILE') AS ref_parent_file;
SELECT getenv('ALT_PARENT_FILE') AS alt_parent_file;
SELECT getenv('MIN_QUAL') AS ref_qual_cutoff;
SELECT getenv('MIN_DEPTH') AS ref_depth_cutoff;
SELECT getenv('REF_OUT_FILE') AS ref_out_file;
SET
    preserve_insertion_order = false;

CREATE
OR REPLACE VIEW parents AS
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
    (
        CASE
            WHEN variant = '.' THEN reference
            ELSE UPPER(variant)
        END
    ) AS mod_variant
FROM
    read_parquet([ref_parent_file, alt_parent_file]);

CREATE
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
    n_variants >= 2;

CREATE
OR REPLACE VIEW temp_ref AS
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
    AND quality > ref_qual_cutoff
    AND LENGTH (mod_variant) <= 1
    AND LENGTH (up_reference) <= 1
    AND depth > ref_depth_cutoff;

CREATE
OR REPLACE VIEW ref_parent AS
SELECT
    sample,
    chromosome,
    int_pos,
    mod_variant AS ref_allele
FROM
    temp_ref
WHERE
    sample = ref_parent_name;

CREATE
OR REPLACE VIEW alt_parent AS
SELECT
    sample,
    chromosome,
    int_pos,
    mod_variant AS alt_allele
FROM
    temp_ref
WHERE
    sample = alt_parent_name;

CREATE
OR REPLACE TABLE ref AS
SELECT
    *
FROM
    ref_parent
    INNER JOIN alt_parent ON (
        ref_parent.chromosome = alt_parent.chromosome
        AND ref_parent.int_pos = alt_parent.int_pos
    );

CREATE
OR REPLACE TABLE refs_unique AS
SELECT
    *
FROM
    ref
WHERE
    ref_allele != alt_allele;

COPY (
    SELECT
        chromosome,
        int_pos AS position,
        ref_allele AS reference,
        alt_allele AS variant
    FROM
        ref
) TO ref_out_file (FORMAT 'parquet');