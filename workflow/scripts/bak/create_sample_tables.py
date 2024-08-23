from pyspark.sql import SparkSession

spark = SparkSession.builder.appName("VCF DB").getOrCreate()

%%sql

CREATE OR REPLACE VIEW temp_vcfs AS
SELECT * FROM read_parquet('{{sample_file_glob}}');

CREATE OR REPLACE VIEW vcfs AS
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
    (CASE WHEN temp_vcfs.variant='.' THEN temp_vcfs.reference ELSE temp_vcfs.variant END) AS new_variant
    FROM temp_vcfs;

CREATE OR REPLACE VIEW samples_rearranged AS
SELECT *, (CASE WHEN ref_allele=new_variant THEN '.' ELSE new_variant END) AS var_adjusted FROM vcfs
    INNER JOIN ref ON vcfs.chromosome = ref.chromosome AND vcfs.int_pos = ref.int_pos
    WHERE new_variant=ref_allele OR new_variant=alt_allele;

CREATE OR REPLACE TABLE final_samples AS
SELECT 
    string_split(sample, '-')[1] AS condition,
    string_split(sample, '-')[2] AS sample_type,
    string_split(sample, '-')[3] AS sample_num,
    reference,
    variant,
    chromosome,
    int_pos AS position,
    string_split(allele_depth, ',')[1] AS ref_reads,
    string_split(allele_depth, ',')[2] AS variant_reads,
    quality AS QUAL,
    genotype AS GT,
    depth AS DP
    FROM samples_rearranged
ORDER BY sample_num, chromosome, position;

CHECKPOINT;