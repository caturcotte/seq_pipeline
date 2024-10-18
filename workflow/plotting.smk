rule format_bcf:
    input:
        multiext("data/calls/{sample}", ".bcf", ".bcf.csi"),
    output:
        temp("data/tsvs/{sample}_tmp.tsv"),
    params:
        format=get_query_format,
    envmodules:
        config["envmodules"]["samtools"],
    shell:
        "bcftools query -f '{params.format}' -o {output} {input[0]}"


rule add_header_to_tsv:
    input:
        "data/tsvs/{sample}_tmp.tsv",
    output:
        "data/tsvs/{sample}.tsv",
    params:
        sedex=lambda w: get_sedex(w),
    shell:
        "sed {params.sedex} {input} > {output}"


rule convert_tsv_to_parquet:
    input:
        "data/tsvs/{sample}.tsv",
    output:
        "data/parquets/{sample}_converted.parquet",
    conda:
        "envs/convert_tsv_to_parquet.yaml"
    params:
        short=lambda w: check_if_short_read(w),
    script:
        "scripts/convert_tsv_to_parquet.py"


rule rearrange_reference_parquet:
    input:
        ref=f"data/parquets/{config['reference_genotype']['name']}_converted.parquet",
        alt=f"data/parquets/{config['variant_genotype']['name']}_converted.parquet",
    output:
        "data/parquets/reference.parquet",
        db=temp("reference.db"),
    params:
        ref_name=config["reference_genotype"]["name"],
        alt_name=config["variant_genotype"]["name"],
        min_qual=config["min_qual"],
        min_depth=config["min_depth"],
    # threads: 16
    shell:
        "ipython workflow/scripts/rearrange_reference_parquet.py -- {input.ref} {input.alt} {params.ref_name} {params.alt_name} {params.min_qual} {params.min_depth} {output[0]} {output.db}"


rule rearrange_progeny_parquet:
    input:
        smp="data/parquets/{progeny}_converted.parquet",
        ref="data/parquets/reference.parquet",
        pat=f"data/parquets/{config['paternal_genotype']['name']}_converted.parquet",
    output:
        "data/parquets/{progeny}_rearranged.parquet",
        temp("{progeny}.db"),
    # threads: 16
    shell:
        "ipython workflow/scripts/rearrange_sample_parquet.py -- {wildcards.progeny} {input.ref} {input.pat} {input.smp} {output[0]}"


rule make_input_output_list_files:
    input:
        expand("data/parquets/{progeny}_rearranged.parquet", progeny=all_progeny),
    output:
        temp("data/parquets/all_samples_input.txt"),
        temp("data/parquets/all_samples_output.txt"),
    run:
        with open(output[0], "w") as file:
            for i in input:
                file.write(f"{i}\n")
        with open(output[1], "w") as file:
            for i in input:
                basename = Path(i).stem.rstrip("_rearranged")
                file.write(f"data/tiger/{basename}_phased.csv\n")


rule exclude_problem_snps:
    input:
        "data/parquets/all_samples_input.txt",
        "data/parquets/all_samples_output.txt",
    output:
        expand("data/plots/{progeny}_phased.csv", progeny=all_progeny),
    params:
        chr_names=list(config["chromosomes"].keys()),
    shell:
        "ipython workflow/scripts/remove_problem_snps.py -- {input[0]} {input[1]} {output[0]} {params.chr_names}"
