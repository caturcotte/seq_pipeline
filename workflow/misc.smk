rule symlink_ref:
    input:
        config["reference"],
    output:
        f"data/resources/{config['ref_name']}.fa",
    shell:
        "ln -s {input} {output}"


rule mask_repeats:
    input:
        "data/resources/{ref}.fa",
    output:
        "data/resources/{ref}.fa.masked",
    threads: 16
    shell:
        "RepeatMasker -species 7227 -pa {threads} -dir resources/ -s {input}"


rule mv_masked_ref:
    input:
        "data/resources/{ref}.fa.masked",
    output:
        "data/resources/{ref}_masked.fa",
    shell:
        "mv {input} {output}"


rule faidx_ref:
    input:
        "{prefix}",
    output:
        "{prefix}.fai",
    # conda:
    #     "envs/samtools.yaml"
    shell:
        "samtools faidx {input}"


rule symlink_fqs_single:
    input:
        get_file_locations,
    output:
        "data/reads/{sample}.fq.gz",
    shell:
        "ln -s {input} {output}"


rule symlink_fqs_paired:
    input:
        r1=lambda w: get_file_locations(w, read="1"),
        r2=lambda w: get_file_locations(w, read="2"),
    output:
        r1="data/reads/{sample}_1.fq.gz",
        r2="data/reads/{sample}_2.fq.gz",
    shell:
        "ln -s {input.r1} {output.r1} && ln -s {input.r2} {output.r2}"


rule trim_adapters:
    input:
        ["data/reads/{sample}_1.fq.gz", "data/reads/{sample}_2.fq.gz"],
    output:
        fastq1="data/reads/{sample}_1_trimmed.fq.gz",
        fastq2="data/reads/{sample}_2_trimmed.fq.gz",
        qc="data/reads/{sample}.qc.txt",
    threads: 2
    params:
        adapters=lambda w: get_adapter_seqs(w),
        extra="--minimum-length 1 -q 20",
    wrapper:
        "v2.6.0/bio/cutadapt/pe"


rule bwa_idx:
    input:
        "{prefix}",
    output:
        multiext("{prefix}", ".ann", ".bwt", ".pac", ".sa", ".amb"),
    conda:
        "envs/bwa.yaml"
    shell:
        "bwa index {input}"


rule minimap2_idx:
    input:
        "{prefix}",
    output:
        "{prefix}.mmi",
    threads: 8
    conda:
        "envs/minimap2.yaml"
    shell:
        "minimap2 -d {input} > {output}"


rule bowtie2_build:
    input:
        "{prefix}",
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
    wrapper:
        "v2.6.0/bio/bowtie2/build"


# rule generate_freebayes_regions:
#     input:
#         ref=get_ref_to_generate_regions,
#         idx=lambda w: get_ref_to_generate_regions(w, fai=True),
#     output:
#         regions=expand("data/resources/regions/{{ref}}.{{chrom}}.region.{i}.bed", i=chunks),
#     params:
#         nchunks=config["freebayes_opts"]["nchunks"],
#     resources:
#         time="5-0",
#     shell:
#         "scripts/fasta_generate_regions.py --chunks --bed=data/resources/regions/{wildcards.ref} --chromosome={wildcards.chrom} {input.idx} {params.nchunks}"


rule bcftools_index:
    input:
        "{prefix}",
    output:
        "{prefix}.csi",
    # conda:
    #     "envs/bcftools.yaml"
    shell:
        "bcftools index {input}"
