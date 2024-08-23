rule symlink_ref:
    input:
        config["reference"],
    output:
        f"data/resources/{config['ref_name']}.fa",
    shell:
        "ln -s {input} {output}"


rule mask_repeats:
    input:
        f"data/resources/{config['ref_name']}.fa",
    output:
        f"data/resources/{config['ref_name']}.fa.masked",
    threads: 16
    cache: True
    envmodules:
        config['envmodules']['repeatmasker']
    shell:
        'RepeatMasker -species "drosophila melanogaster" -pa {threads} -dir data/resources/ -s {input}'


rule mv_masked_ref:
    input:
        f"data/resources/{config['ref_name']}.fa.masked",
    output:
        f"data/resources/{config['ref_name']}_masked.fa",
    shell:
        "mv {input} {output}"


rule faidx_ref:
    input:
        "{prefix}",
    output:
        "{prefix}.fai",
    cache: True
    envmodules:
        config['envmodules']['samtools']
    shell:
        "samtools faidx {input}"


rule symlink_fastqs:
    input:
        get_fastq_from_id,
    output:
        "data/reads/{sample}_{iden}{read}.fq.gz",
    shell:
        "ln -s {input} {output}"


rule minimap2_idx:
    input:
        "{prefix}",
    output:
        "{prefix}.mmi",
    threads: 8
    envmodules:
        config['envmodules']['minimap2']
    cache: True
    shell:
        "minimap2 -ax map-ont -d {output} {input}"


rule bowtie2_build:
    input:
        ref=get_ref,
    output:
        multiext(
            "{prefix}",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    cache: True
    envmodules:
        config['envmodules']['bowtie2']
    params:
        base=lambda w: get_ref(w.prefix, base=True)
    shell:
        "bowtie2-build {input.ref} {params.base}"


rule bcftools_index:
    input:
        "{prefix}",
    output:
        "{prefix}.csi",
    envmodules:
        config['envmodules']['samtools']
    shell:
        "bcftools index {input}"


rule gzip_fastq:
    input:
        "{prefix}.fastq",
    output:
        "{prefix}.fastq.gz"
    shell:
        "gzip {input}"