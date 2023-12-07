rule bcftools_mpileup_single:
    input:
        alns="data/alignments/{sample}_dedup.bam",
        idxs="data/alignments/{sample}_dedup.bam.bai",
        ref=get_ref,
        ref_idx=lambda w: get_ref(w, fai=True),
    output:
        bcf=temp("data/calls/{sample}_bcftools_pileup.bcf"),
        call_type=temp(".tmp/single_call_{sample}.txt"),
    params:
        extra="--max-depth 200 --min-BQ 15",
    threads: 16
    shell:
        "bcftools mpileup -f {input.ref} {params.extra} -Ou -o {output.bcf} --threads {threads} {input.alns} && touch {output.call_type}"


rule bcftools_mpileup_group:
    input:
        alns=get_alns_in_group,
        idxs=lambda w: get_alns_in_group(w, bai=True),
        ref=get_ref,
        ref_idx=lambda w: get_ref(w, fai=True),
    output:
        "data/calls/{group}_bcftools_pileup.bcf",
    params:
        extra="--max-depth 200 --min-BQ 15",
    threads: 16
    shell:
        "bcftools mpileup -f {input.ref} {params.extra} -Ou -o {output} --threads {threads} {input.alns}"


rule bcftools_call:
    input:
        "data/calls/{sample_or_group}_bcftools_pileup.bcf",
    output:
        bcf="data/calls/{sample_or_group}_bcftools_unprocessed.bcf",
    threads: 4
    shell:
        "bcftools call -cv -o {output} --threads {threads} {input}"


ruleorder: bcftools_call > separate_into_samples


rule freebayes:
    input:
        alns=get_alns_in_group,
        idxs=lambda w: get_alns_in_group(w, bai=True),
        region=get_region_from_sample,
        ref=get_ref,
        ref_idx=lambda w: get_ref(w, fai=True),
    output:
        "data/calls/{group}_{label}_freebayes_{i}.bcf",
    threads: 1
    resources:
        time="2-0",
        mem_mb="1000",
    shell:
        "freebayes -f {input.ref} -t {input.region} {input.alns} | bcftools view -Ob -o {output}"

rule concat_freebayes:
    input:
        calls=expand("data/calls/{{group}}_{label}_freebayes_{i}.bcf", i=chunks, label=get_labels()),
        idxs=expand("data/calls/{{group}}_{label}_freebayes_{i}.bcf.csi", i=chunks, label=get_labels()),
    output:
        "data/calls/{group}_freebayes_unprocessed.bcf",
    threads: 4
    shell:
        "bcftools concat -a -D {input.calls} -o {output}"


rule separate_into_samples:
    input:
        get_group_from_sample,
    output:
        bcfs="data/calls/{sample}_{caller}_unprocessed.bcf",
        call_type=temp(".tmp/group_call_{sample}_{caller}.txt"),
    resources:
        time="2:00:00",
    shell:
        "bcftools view -s {wildcards.sample} -a -o {output.bcfs} {input} && touch {output.call_type}"
