rule convert_tsv_to_parquet:
    input:
        "data/tsvs/{sample}.tsv",
    output:
        "data/parquets/{sample}_converted.parquet",
    conda:
        "envs/convert_tsv_to_parquet.yaml"
    script:
        "scripts/convert_tsv_to_parquet.py"


rule rearrange_reference_parquet:
    input:
        ref=f"data/parquets/{config['ref_parent']}_converted.parquet",
        alt=f"data/parquets/{config['alt_parent']}_converted.parquet",
    output:
        "data/parquets/reference.parquet",
        db=temp("data/resources/default.db"),
    params:
        min_qual=config["min_qual"],
        min_depth=config["min_depth"],
    # conda:
    #     "envs/rearrange_parquet.yaml"
    threads: 16
    shell:
        """
        export REF_PARENT_FILE={input.ref} &&
        export ALT_PARENT_FILE={input.alt} &&
        export MIN_QUAL={params.min_qual} &&
        export MIN_DEPTH={params.min_depth} &&
        export REF_OUT_FILE={output[0]} &&
        workflow/scripts/duckdb {output.db} < workflow/scripts/rearrange_reference_parquet.sql
        """


rule rearrange_sample_parquet:
    input:
        smp="data/parquets/{sample}_converted.parquet",
        ref="data/parquets/reference.parquet",
    output:
        "data/parquets/{sample}.parquet",
        db=temp("data/resources/{sample}.db"),
    conda:
        "envs/rearrange_parquet.yaml"
    threads: 16
    script:
        "scripts/rearrange_sample_parquet.py"


rule make_tiger_input:
    input:
        "data/parquets/{sample}.parquet",
    output:
        "data/tiger/{sample}_tiger_input.txt",
    params:
        chr_names=config["chr_names"],
    run:
        tiger_chrom_name_dict = dict(
            zip(params["chr_names"], range(len(params["chr_names"])))
        )
        df = read_parquet(input[0])
        df["tiger_chrom"] = df.apply(
            lambda row: tiger_chrom_name_dict[row["chromosome"]], axis=1
        )
        df.to_csv(
            output[0],
            sep="\t",
            header=False,
            index=False,
            columns=[
                "tiger_chrom",
                "position",
                "reference",
                "ref_reads",
                "variant",
                "variant_reads",
            ],
        )


rule run_tiger:
    input:
        "data/tiger/{sample}_tiger_input.txt",
    output:
        "data/tiger/{sample}/{sample}_CO_estimates.txt",
    envmodules:
        config["envmodules"]["r"],
        config["envmodules"]["perl"],
    shell:
        "sh scripts/run_tiger.sh {input} data/tiger/{wildcards.sample}"


rule process_tiger_output:
    input:
        expand("data/tiger/{sample}/{sample}_CO_estimates.txt", sample=progeny),
    output:
        marker_df="data/tiger/marker_df.pq",
        pre_intervals_df="data/tiger/pre_intervals_df.pq",
        intervals_df="data/tiger/intervals_df.pq",
    conda:
        "envs/process_tiger_outputs.yaml"
    script:
        "scripts/process_tiger_output.py"
