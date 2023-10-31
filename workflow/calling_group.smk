rule bcftools_group_call:
    input:
        alns=get_alns_for_calling,
        aln_idxs=lambda w: get_alns_for_calling(w, bai=True),
        ref=get_ref,
        ref_idx=lambda w: get_ref(w, fai=True),
    output:
        "called/{group}/bcftools_unprocessed.bcf",
    shell:
        "bcftools call -cv -f {input.ref} {input.alns} -o {output}"


rule freebayes:
    input:
        alns=get_alns_for_calling,
        idxs=lambda w: get_alns_for_calling(w, bai=True),
        region=get_region_from_sample,
        ref=get_ref,
        ref_idx=lambda w: get_ref(w, fai=True),
    output:
        temp("called/{chrom}/{group}_{i}.bcf"),
        temp("called/{chrom}/{group}_{i}.bcf.csi"),
    log:
        "logs/freebayes/{chrom}_{group}_{i}.log",
    threads: 1
    resources:
        mem_mb=50000,
        time="1-0",
    shell:
        "freebayes -f {input.ref} -t {input.region} {input.alns} | bcftools view -Ob -o {output} --write-index"


rule concat_freebayes:
    input:
        calls=expand("called/{chrom}/{{group}}_{i}.bcf", i=chunks, chrom=chroms),
        idxs=expand("called/{chrom}/{{group}}_{i}.bcf.csi", i=chunks, chrom=chroms),
    output:
        temp("called/{group}/freebayes_unprocessed.bcf"),
    log:
        "logs/concat_vcfs/{group}.log",
    threads: 4
    conda:
        "envs/concat_vcfs.yaml"
    shell:
        "bcftools concat -a {input.calls} | vcfuniq > {output} 2> {log}"


rule separate_into_samples:
    input:
        get_group_from_sample,
    output:
        "called/{sample}_{caller}_unprocessed.bcf",
    resources:
        time="2:00:00",
    shell:
        "bcftools view -s {wildcards.sample} -a -o {output} {input}"
