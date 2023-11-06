rule symlink_ref:
    input:
        config["reference"],
    output:
        f"resources/{config['ref_name']}.fa",
    shell:
        "ln -s {input} {output}"


rule mask_repeats:
    input:
        "resources/{ref}.fa",
    output:
        "resources/{ref}.fa.masked",
    threads: 16
    shell:
        "RepeatMasker -species 7227 -pa {threads} -dir resources/ -s {input}"


rule faidx_ref:
    input:
        "{prefix}",
    output:
        "{prefix}.fai",
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools faidx {input}"


rule symlink_fqs_single:
    input:
        get_file_locations,
    output:
        "reads/{sample}.fq.gz",
    shell:
        "ln -s {input} {output}"


rule symlink_fqs_paired:
    input:
        r1=lambda w: get_file_locations(w, end="r1"),
        r2=lambda w: get_file_locations(w, end="r2"),
    output:
        r1="reads/{sample}_1.fq.gz",
        r2="reads/{sample}_2.fq.gz",
    shell:
        "ln -s {input.r1} {output.r1} && ln -s {input.r2} {output.r2}"


rule trim_adapters:
    input:
        ["reads/{sample}_1.fq.gz", "reads/{sample}_2.fq.gz"],
    output:
        fastq1="reads/{sample}_1_unsorted.fq.gz",
        fastq2="reads/{sample}_2_unsorted.fq.gz",
        qc="reads/{sample}.qc.txt",
    threads: 8
    params:
        adapters=f"-a {config['trimming_opts']['adapters']['fwd']} -g {config['trimming_opts']['adapters']['rev']} -A {config['trimming_opts']['adapters']['fwd']} -G {config['trimming_opts']['adapters']['rev']}",
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


rule generate_freebayes_regions:
    input:
        ref=get_ref_to_generate_regions,
        idx=lambda w: get_ref_to_generate_regions(w, fai=True),
    output:
        regions=expand("resources/regions/{{ref}}.{{chrom}}.region.{i}.bed", i=chunks),
    params:
        nchunks=config["freebayes_opts"]["nchunks"],
    resources:
        time="5-0",
    shell:
        "scripts/fasta_generate_regions.py --chunks --bed=resources/regions/{wildcards.ref} --chromosome={wildcards.chrom} {input.idx} {params.nchunks}"


rule bcftools_index:
    input:
        "called/{sample}_{caller}_unprocessed.bcf",
    output:
        "called/{sample}_{caller}_unprocessed.bcf.csi",
    conda:
        "envs/bcftools.yaml"
    shell:
        "bcftools index {input}"


rule vcf_stats:
    input:
        bcf="called/{sample}.bcf",
        ref=get_ref,
    output:
        "metrics/{sample}.bcf.stats",
    conda:
        "envs/bcftools.yaml"
    shell:
        "bcftools stats -F {input.ref} -s - {input.bcf} > {output}"
