def get_file_locations(w, end=None):
    sample = sample_sheet.loc[sample_sheet["sample"] == w.sample]
    if config["paired"]:
        if end == "r1":
            return sample["path_to_reads_1"]
        elif end == "r2":
            return sample["path_to_reads_2"]
    else:
        return sample["path_to_reads"]


rule symlink_ref:
    input:
        config["reference"],
    output:
        f"resources/{config['ref_name']}.fa",
    shell:
        "ln -s {input} {output}"

rule symlink_ref_idx:
    input:
        config["reference"] + ".fai",
    output:
        f"resources/{config['ref_name']}.fa.fai",
    shell:
        "ln -s {input} {output}"

rule symlink_fqs_single:
    input:
        get_file_locations,
    output:
        "reads/{sample}.fq.gz",
    shell:
        "ln -s {input} {output}"

rule generate_freebayes_regions:
    input:
        ref=get_ref_to_generate_regions,
        idx=lambda w: get_ref_to_generate_regions(w, fai=True)
    output:
        regions=expand("resources/regions/{{ref}}.{{chrom}}.region.{i}.bed", i=chunks),
    params:
        nchunks=config["freebayes_opts"]["nchunks"],
    resources:
        time="5-0",
        mem_mb=100000,
    shell:
        "scripts/fasta_generate_regions.py --chunks --bed=resources/regions/{wildcards.ref} --chromosome={wildcards.chrom} {input.idx} {params.nchunks}"
rule bcftools_index:
    input:
        "called/{sample}_{caller}_unprocessed.bcf",
    output:
        "called/{sample}_{caller}_unprocessed.bcf.csi",
    shell:
        "bcftools index {input}"


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
    resources:
        mem_mb=10000,
    params:
        adapters=f"-a {config['trimming_opts']['adapters']['fwd']} -g {config['trimming_opts']['adapters']['rev']} -A {config['trimming_opts']['adapters']['fwd']} -G {config['trimming_opts']['adapters']['rev']}",
        extra="--minimum-length 1 -q 20",
    wrapper:
        "v2.3.1/bio/cutadapt/pe"


rule mask_repeats:
    input:
        "resources/{ref}.fa",
    output:
        "resources/{ref}.fa.masked",
    threads: 8
    shell:
        "RepeatMasker -species drosophila -pa {threads} -dir resources/ -s {input}"

rule minimap2_idx:
    input:
        "{prefix}",
    output:
        "{prefix}.mmi"
    threads: 8
    shell:
        "minimap2 -d {input} > {output}"

rule bowtie2_build:
    input:
        "resources/{ref}.fa.masked",
    output:
        multiext(
            "resources/{ref}.fa.masked",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    threads: 8
    wrapper:
        "v2.6.0/bio/bowtie2/build"


rule faidx_ref:
    input:
        "resources/{ref}.fa.masked",
    output:
        "resources/{ref}.fa.masked.fai",
    shell:
        "samtools faidx {input}"


rule vcf_stats:
    input:
        bcf="called/{sample}.bcf",
        ref=reference,
    output:
        "metrics/{sample}.bcf.stats",
    shell:
        "bcftools stats -F {input.ref} -s - {input.bcf} > {output}"
