# for some reason you can't just select for snps and the reference so instead you just have to exclude everything else
# this weirdly doesn't get rid of all indels so I remove them when rearranging the parquets
rule select_snps:
    input:
        bcf="data/calls/{sample}_raw.bcf",
        csi="data/calls/{sample}_raw.bcf.csi",
    output:
        "data/calls/{sample}.bcf",
    threads: 4
    envmodules:
        config['envmodules']['samtools']
    shell:
        "bcftools view -V mnps,bnd,indels,other -Ob -o {output} {input.bcf}"


rule format_bcf:
    input:
        bcf="data/calls/{sample}.bcf",
        csi="data/calls/{sample}.bcf.csi",
    output:
        temp("data/tsvs/{sample}_tmp.tsv"),
    params:
        format=get_query_format,
    envmodules:
        config['envmodules']['samtools']
    shell:
        "bcftools query -f '{params.format}' -o {output} {input.bcf}"


rule add_header_to_tsv:
    input:
        "data/tsvs/{sample}_tmp.tsv",
    output:
        "data/tsvs/{sample}.tsv",
    params:
        sedex="\'1s/^/sample\\tchromosome\\tposition\\treference\\tvariant\\tquality\\tgenotype\\tdepth\\tallele_depth\\n/\'"
    shell:
        "sed {params.sedex} {input} > {output}"


rule convert_tsv_to_parquet:
    input:
        "data/tsvs/{sample}.tsv",
    output:
        "data/parquets/{sample}.parquet",
    conda:
        "envs/convert_tsv_to_parquet.yaml"
    resources:
        mem_mb=50000,
    script:
        "scripts/convert_tsv_to_parquet.py"

# rearrange the tables so that the reference is ref_parent
rule rearrange_parquet:
    input:
        ref=f"data/parquets/{config['ref_parent']}.parquet",
        alt=f"data/parquets/{config['alt_parent']}.parquet",
        prog=expand("data/parquets/{sample}.parquet", sample=progeny),
    output:
        ref="data/parquets/reference.parquet",
        prog="data/parquets/samples.parquet",
    conda:
        "envs/rearrange_parquet.yaml"
    notebook:
        "scripts/rearrange_parquet.py.ipynb"


rule run_tiger:
    input:
        "data/parquets/reference.parquet",
        "data/parquets/samples.parquet",
    output:
        "data/tiger/tiger_output.parquet",
    conda:
        "envs/run_tiger.yaml"
    script:
        "scripts/run_tiger.py"
