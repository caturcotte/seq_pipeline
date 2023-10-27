rule bcftools_call_single:
    input:
        alns="mapped/{sample}_sort_dedup.bam",
        aln_idxs="mapped/{sample}_sort_dedup.bam.bai",
        ref=get_ref,
        ref_idx=lambda w: get_ref(w, fai=True),
    output:
        "called/{sample}_bcftools_unprocessed.bcf",
    shell:
        "bcftools call -cv -f {input.ref} {input.alns} -o {output}"


rule freebayes_single:
    input:
        alns="mapped/{sample}_sort_dedup.bam",
        aln_idxs="mapped/{sample}_sort_dedup.bam.bai",
        region=get_region_from_sample,
        ref=get_ref,
        ref_idx=lambda w: get_ref(w, fai=True),
    output:
        temp("called/{chrom}/{sample}_{i}.bcf"),
        temp("called/{chrom}/{sample}_{i}.bcf.csi"),
    log:
        "logs/freebayes/{chrom}_{sample}_{i}.log",
    threads: 1
    resources:
        mem_mb=50000,
        time="1-0",
    shell:
        "freebayes -f {input.ref} -t {input.region} {input.alns} | bcftools view -Ob -o {output} --write-index"

rule concat_freebayes_vcfs_single:
    input:
        calls=expand("called/{chrom}/{{sample}}_{i}.bcf", i=chunks, chrom=chroms),
        idxs=expand("called/{chrom}/{{sample}}_{i}.bcf.csi", i=chunks, chrom=chroms)
    output:
        temp("called/{sample}_freebayes_unprocessed.bcf"),
    log:
        "logs/concat_vcfs/{sample}.log",
    threads: 4
    conda:
        "envs/concat_vcfs.yaml"
    shell:
        "bcftools concat -a {input.calls} | vcfuniq > {output} 2> {log}"

