rule bcftools_mpileup_single:
    input:
        alns="data/alignments/{sample}_dedup.bam",
        idxs="data/alignments/{sample}_dedup.bam.bai",
        ref=get_ref,
        ref_idx=lambda w: get_ref(w, fai=True),
    output:
        bcf=temp("data/calls/{sample}_bcftools_single_pileup.bcf"),
    params:
        extra=config["bcftools_opts"]["pileup_options"],
    threads: 16
    shell:
        "bcftools mpileup -f {input.ref} -Ou -o {output.bcf} {params.extra} --threads {threads} {input.alns}"


rule bcftools_mpileup_group:
    input:
        alns=get_alns_in_group,
        idxs=lambda w: get_alns_in_group(w, bai=True),
        ref=get_ref,
        ref_idx=lambda w: get_ref(w, fai=True),
    output:
        bcf="data/calls/{group}_bcftools_group_pileup.bcf",
    params:
        extra=config["bcftools_opts"]["pileup_options"],
    resources:
        mem_mb=50000,
        time="2-0",
    threads: 32
    shell:
        "bcftools mpileup -f {input.ref} -Ou -o {output.bcf} {params.extra} --threads {threads} {input.alns}"


rule bcftools_call:
    input:
        pileup=lambda w: group_or_single(w),
    output:
        bcf="data/calls/{sample_or_group}_bcftools_unprocessed.bcf",
        txt=temp(".tmp/{sample_or_group}_bcftools.txt"),
    threads: 16
    params:
        call_type=config["bcftools_opts"]["call_type"],
    resources:
        time="1-0",
    shell:
        "bcftools call -{params.call_type}v -o {output.bcf} --threads {threads} {input.pileup} && touch {output.txt}"


rule freebayes:
    input:
        alns=get_alns_in_group,
        idxs=lambda w: get_alns_in_group(w, bai=True),
        region=get_region_from_sample,
        ref=get_ref,
        ref_idx=lambda w: get_ref(w, fai=True),
    output:
        bcf=temp("data/calls/{group}_{label}_freebayes_{i}.bcf"),
    threads: 1
    resources:
        time="2-0",
        mem_mb="1000",
    shell:
        "freebayes -f {input.ref} -t {input.region} {input.alns} | bcftools view -Ob -o {output.bcf}"


rule concat_freebayes:
    input:
        calls=expand(
            "data/calls/{{group}}_{label}_freebayes_{i}.bcf",
            i=chunks,
            label=get_labels(),
        ),
        idxs=expand(
            "data/calls/{{group}}_{label}_freebayes_{i}.bcf.csi",
            i=chunks,
            label=get_labels(),
        ),
    output:
        bcf="data/calls/{group}_freebayes_unprocessed.bcf",
        txt=temp(".tmp/group_call_{group}_freebayes.txt"),
    threads: 4
    shell:
        "bcftools concat -a -D {input.calls} -o {output.bcf} && touch {output.txt}"


rule separate_into_samples:
    input:
        get_group_from_sample,
    output:
        bcfs="data/calls/{sample}_{caller}_unprocessed.bcf",
    resources:
        time="2:00:00",
    shell:
        "bcftools view -s {wildcards.sample} -a -o {output.bcfs} {input}"


ruleorder: bcftools_call > separate_into_samples
