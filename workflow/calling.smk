rule bcftools_mpileup_single:
    input:
        alns="data/alignments/{sample}_markdup.bam",
        idxs="data/alignments/{sample}_markdup.bam.bai",
        ref=get_ref,
        ref_idx=lambda w: get_ref(w, fai=True),
    output:
        bcf=temp("data/calls/{sample}_bcftools_single_pileup.bcf"),
    params:
        extra=config["calling"]["bcftools_opts"]["pileup_options"],
        # regions=lambda w: get_regions(w),
    resources:
        mem_mb=50000,
        time="2-0",
    threads: 32
    shell:
        "bcftools mpileup -f {input.ref} -a AD,DP,SP -Ou -o {output.bcf} {params.extra} --threads {threads} {input.alns}"


rule bcftools_mpileup_group:
    input:
        alns=get_alns_in_group,
        idxs=lambda w: get_alns_in_group(w, bai=True),
        ref=get_ref,
        ref_idx=lambda w: get_ref(w, fai=True),
    output:
        bcf=temp("data/calls/{group}_bcftools_group_pileup.bcf"),
    params:
        extra=config["calling"]["bcftools_opts"]["pileup_options"],
        # regions=get_regions(),
    resources:
        mem_mb=50000,
        time="2-0",
    threads: 128
    shell:
        "bcftools mpileup -f {input.ref} -a AD,DP,SP -Ou -o {output.bcf} {params.extra} --threads {threads} {input.alns}"


rule bcftools_call:
    input:
        pileup=lambda w: group_or_single(w),
    output:
        bcf="data/calls/{sample_or_group}_bcftools_raw.bcf",
        txt=temp(".tmp/{sample_or_group}_bcftools.txt"),
    threads: 128
    params:
        call_type=config["calling"]["bcftools_opts"]["call_type"],
        # regions=get_regions(),
    resources:
        mem_mb=50000,
        time="2-0",
    shell:
        "bcftools call -{params.call_type}v -f GQ,GP -o {output.bcf} --threads {threads} {input.pileup} && touch {output.txt}"


rule freebayes:
    input:
        alns=get_alns_in_group,
        idxs=lambda w: get_alns_in_group(w, bai=True),
        region=get_region_from_sample,
        ref=get_ref,
        ref_idx=lambda w: get_ref(w, fai=True),
    output:
        bcf=temp("data/calls/{group}_{label}_freebayes_{i}.bcf"),
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
        bcf="data/calls/{group}_freebayes_raw.bcf",
        txt=temp(".tmp/group_call_{group}_freebayes.txt"),
    threads: 4
    shell:
        "bcftools concat -a -D {input.calls} -o {output.bcf} && touch {output.txt}"


rule overlap_bcftools_and_freebayes:
    input:
        fb="data/calls/{group}_freebayes_raw.bcf",
        fb_idx="data/calls/{group}_freebayes_raw.bcf.csi",
        bt="data/calls/{group}_bcftools_raw.bcf",
        bt_idx="data/calls/{group}_bcftools_raw.bcf.csi",
    output:
        temp(expand("data/calls/bt_fb_overlap_{{group}}/000{i}.vcf", i=[str(i) for i in range(4)])),
    shell:
        "bcftools isec -n 2 -p bt_fb_overlap_{wildcards.group} {input.fb} {input.bt}"

rule annotate_freebayes_with_bcftools:
    input:
        fb="data/calls/bt_fb_overlap_{group}/0000.vcf",
        bt="data/calls/{group}_bcftools_raw.bcf",
    output:
        "data/calls/{group}_both_raw.bcf",
    shell:
        "bcftools annotate -a {input.bt} -c 'CHROM,POS,REF,ALT,=INFO,=FORMAT' -l unique -Ob -o {output} {input.fb}"
