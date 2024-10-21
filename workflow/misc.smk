from pathlib import Path


rule download_dorado:
    output:
        f"workflow/tools/dorado/dorado-0.8.1-{os}-{arch}/bin/dorado",
    params:
        out_dir=lambda w, output: str(Path(output[0]).parents[2]),
    shell:
        "wget -qO- https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.8.1-{os}-{arch}.tar.gz | tar xvz -C {params.out_dir}"


rule download_rerio:
    output:
        "data/resources/rerio/download_model.py",
    params:
        out_dir=lambda w, output: str(Path(output[0]).parents[0]),
    shell:
        "git clone https://github.com/nanoporetech/rerio {params.out_dir}"


rule download_clair3_model:
    input:
        "data/resources/rerio/download_model.py",
    output:
        "data/resources/rerio/clair3_models/r1041_e82_400bps_sup_v500/pileup.index",
    params:
        base_dir=lambda w, input: str(Path(input[0]).parents[0].joinpath('clair3_models/r1041_e82_400bps_sup_v500_model'))
    shell:
        "{input} {params.base_dir}"


rule symlink_ref:
    input:
        config["reference"]['location'],
    output:
        f"data/resources/genome.fa",
    shell:
        "ln -s {input} {output}"


rule index_bam:
    input:
        "{prefix}.bam",
    output:
        "{prefix}.bam.bai",
    envmodules:
        config['envmodules']['sambamba']
    shell:
        "sambamba index {input}"


# rule vcf_to_bcf:
#     input:
#         "{prefix}.vcf",
#     output:
#         "{prefix}.bcf",
#     shell:
#         "bcftools view -Ob -o {output} {input}"


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
        "{prefix}.fa",
    output:
        "{prefix}.mmi",
    threads: 8
    envmodules:
        config['envmodules']['minimap2']
    # cache: True
    shell:
        "minimap2 -x lr:hq -d {output} {input}"


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
        base=lambda w: get_ref(w.prefix, base=True),
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
        "{prefix}.fastq.gz",
    shell:
        "gzip {input}"
