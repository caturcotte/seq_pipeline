rule bcftools_mpileup_single:
    input:
        alns='mapped/{sample}_sort_dedup.bam',
        idxs='mapped/{sample}_sort_dedup.bam.bai',
        ref=get_ref,
        ref_idx=lambda w: get_ref(w, fai=True),
    output:
        bcf=temp("called/{sample}_bcftools_pileup.bcf"),
        call_type=temp('.tmp/single_call_{sample}.txt')
    params:
        regions=lambda w: get_regions_to_call(w),
        extra="--max-depth 200 --min-BQ 15"
    threads: 16
    # conda:
    #     "envs/bcftools.yaml"
    shell:
        "bcftools mpileup -f {input.ref} -r {params.regions} {params.extra} -Ou -o {output.bcf} --threads {threads} {input.alns} && touch {output.call_type}"


rule bcftools_mpileup_group:
    input:
        alns=get_alns_in_group,
        idxs=lambda w: get_alns_in_group(w, bai=True),
        ref=get_ref,
        ref_idx=lambda w: get_ref(w, fai=True),
    output:
        "called/{group}_bcftools_pileup.bcf",
    params:
        regions=lambda w: get_regions_to_call(w),
        extra="--max-depth 200 --min-BQ 15"
    threads: 16
    # conda:
    #     "envs/bcftools.yaml"
    shell:
        "bcftools mpileup -f {input.ref} -r {params.regions} {params.extra} -Ou -o {output} --threads {threads} {input.alns}"


rule bcftools_call:
    input:
        "called/{sample_or_group}_bcftools_pileup.bcf",
    output:
        bcf="called/{sample_or_group}_bcftools_unprocessed.bcf",
    params:
        # regions=lambda w: get_regions_to_call(w),
    threads: 4
    # conda:
        # "envs/bcftools.yaml"
    shell:
        "bcftools call -cv -o {output} --threads {threads} {input}"


rule freebayes:
    input:
        alns=get_alns_for_pileup,
        idxs=lambda w: get_alns_for_pileup(w, bai=True),
        region=get_region_from_sample,
        ref=get_ref,
        ref_idx=lambda w: get_ref(w, fai=True),
    output:
        temp("called/{chrom}/{group}_{i}.bcf"),
    log:
        "logs/freebayes/{chrom}_{group}_{i}.log",
    threads: 1
    resources:
        time="1-0",
    conda:
        "envs/freebayes.yaml"
    shell:
        "freebayes -f {input.ref} -t {input.region} {input.alns} | bcftools view -Ob -o {output}"


rule concat_freebayes:
    input:
        calls=expand("called/{chrom}/{{group}}_{i}.bcf", i=chunks, chrom=chroms),
        idxs=expand("called/{chrom}/{{group}}_{i}.bcf.csi", i=chunks, chrom=chroms),
    output:
        temp("called/{group}/freebayes_unprocessed.bcf"),
    log:
        "logs/concat_vcfs/{group}.log",
    conda:
        "envs/concat_vcfs.yaml"
    shell:
        "bcftools concat -a {input.calls} | vcfuniq > {output} 2> {log}"


rule separate_into_samples:
    input:
        get_group_from_sample,
    output:
        bcfs="called/{sample}_{caller}_unprocessed.bcf",
        call_type=temp('.tmp/group_call_{sample}_{caller}.txt')
    resources:
        time="2:00:00",
    # conda:
    #     "envs/bcftools.yaml"
    shell:
        "bcftools view -s {wildcards.sample} -a -o {output} {input} && touch {output.call_type}"
