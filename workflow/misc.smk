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
    shell:
        "RepeatMasker -species 7227 -pa {threads} -dir data/resources/ -s {input}"


# rule write_regions:
#     output:
#         temp("data/resources/regions.txt"),
#     params:
#         chromosomes=config["calling"]["chromosomes"],
#         regions=config["calling"]["regions"],
#     run:
#         with open({output}, "w") as file:
#             for i in {params.chromosomes}:
#                 file.write(f"{i}:0-0")
#             for i in {params.regions}:
#                 file.write(i)


# rule create_bed_file:
#     input:
#         "data/resources/regions.txt",
#     output:
#         "data/resources/call_regions.bed",
#     params:
#         regions=config["calling"]["regions"],
#     shell:
#         "sed -e 's/:/\t/' -e 's/-/\t' -i {input} > {output}"


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
    shell:
        "samtools faidx {input}"


rule symlink_fqs_single:
    input:
        get_reads,
    output:
        "data/reads/{sample}_{iden}.fq.gz",
    shell:
        "ln -s {input} {output}"


rule symlink_fqs_paired:
    input:
        get_reads,
    output:
        ["data/reads/{sample}_{iden}_1.fq.gz", "data/reads/{sample}_{iden}_2.fq.gz"],
    shell:
        "ln -s {input[0]} {output[0]} && ln -s {input[1]} {output[1]}"


rule trim_adapters_paired:
    input:
        ["data/reads/{sample}_{prefix}_1.fq.gz", "data/reads/{sample}_{prefix}_2.fq.gz"],
    output:
        fastq1="data/reads/{sample}_{prefix}_1_trimmed.fq.gz",
        fastq2="data/reads/{sample}_{prefix}_2_trimmed.fq.gz",
        qc="data/qc/cutadapt/{sample}_{prefix}.txt",
    threads: 2
    params:
        adapters=lambda w: get_adapter_seqs(w),
        extra="--minimum-length 1 -q 20",
    wrapper:
        "v2.6.0/bio/cutadapt/pe"


rule trim_adapters_single:
    input:
        "data/reads/{sample}_{prefix}.fq.gz",
    output:
        fastq="data/reads/{sample}_{prefix}_trimmed.fq.gz",
        qc="data/qc/cutadapt/{sample}_{prefix}.txt",
    threads: 2
    params:
        adapters=lambda w: get_adapter_seqs(w),
        extra="--minimum-length 1 -q 20",
    wrapper:
        "v2.6.0/bio/cutadapt/se"


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
        ref="{prefix}",
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


# rule generate_chroms_bed:
#     input:
#         get_ref(fai=True),
#     output:
#         "data/resources/call_chromosomes.bed"
#     params:
#         chroms=config['calling']['chromosomes']
#     shell:
#         "awk 'BEGIN {FS=\"\t\"}; {print $1 FS \"0\" FS $2}' {input} > {output}"
# rule merge_chrom_and_region_beds:
#     input:
#         "data/resources/call_chromosomes.bed",
#         "data/resources/call_regions.bed",
#     output:
#         "data/resources/call_chromosomes_and_regions.bed",
#     shell:
#         "cat {input} > {output}"


rule get_fasta_for_calling_regions:
    input:
        bed="data/resources/call_regions.bed",
        ref=get_ref,
    output:
        "data/resources/call_regions.fa",
    shell:
        "bedtools getfasta -fi {input.ref} -bed {input.bed} -fo {output}"


rule generate_freebayes_regions:
    input:
        ref=get_ref,
        idx=lambda w: get_ref(w, fai=True),
    output:
        regions=expand(
            f"data/resources/regions/{config['ref_name']}." + "{label}.region.{i}.bed",
            i=chunks,
            label=get_labels(),
        ),
    params:
        nchunks=config["calling"]["freebayes_opts"]["chunks"],
        basename=config["ref_name"],
    resources:
        time="5-0",
    shell:
        "scripts/fasta_generate_regions.py --chunks --bed=data/resources/regions/{params.basename} {input.idx} {params.nchunks}"


rule bcftools_index:
    input:
        "{prefix}",
    output:
        "{prefix}.csi",
    # conda:
    #     "envs/bcftools.yaml"
    shell:
        "bcftools index {input}"


rule bcf_to_vcf_gz:
    input:
        "{prefix}.bcf",
    output:
        "{prefix}.vcf.gz",
    shell:
        "bcftools view -Oz {input} {output}"


rule separate_ref_fasta:
    input:
        f"data/resources/{config['ref_name']}.fa",
    output:
        expand(
            f"data/resources/{config['ref_name']}.fa.split/{config['ref_name']}"
            + ".part_{label}.fa",
            label=get_labels(),
        ),
    shell:
        "seqkit split -i --two-pass --force {input}"


rule mv_fastas:
    input:
        get_fasta_to_mv,
    output:
        "data/consensus/{sample_or_ref}_{label}.fa",
    shell:
        "mv {input} {output}"


rule separate_sample_fastas:
    input:
        "data/consensus/{sample}.fa",
    output:
        expand(
            "data/consensus/{{sample}}.fa.split/{{sample}}.part_{label}.fasta",
            label=get_labels(),
        ),
    shell:
        "seqkit split -i --two-pass --force {input}"


rule edit_fasta_id:
    input:
        "data/consensus/{sample_or_ref}_{label}.fa",
    output:
        "data/consensus/{sample_or_ref}_{label}_idfix.fa",
    shell:
        "sed -r 's/^(>)/\\1{wildcards.sample_or_ref}_/' {input} > {output}"


rule cat_seqs:
    input:
        get_fas_to_cat,
    output:
        "data/consensus/{group}_{label}.fa",
    shell:
        "cat {input} > {output}"


rule multiple_alignment:
    input:
        "data/consensus/{group}_{label}.fa",
    output:
        "data/msa/{label}_{group}.afa",
    threads: 4
    shell:
        "muscle -align {input} -output {output} -threads {threads}"
