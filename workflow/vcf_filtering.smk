rule bcftools_norm:
    input:
        calls=get_caller,
        ref=get_ref,
        call_type=get_call_type,
    output:
        bcf=temp("called/{sample}_norm.bcf"),
    resources:
        time="2:00:00",
    # conda:
    #     "envs/bcftools.yaml"
    shell:
        "bcftools norm -f {input.ref} -m '-' -Ob -o {output.bcf} {input.calls}"


rule filter_low_quality:
    input:
        bcf="called/{sample}_norm.bcf",
        csi="called/{sample}_norm.bcf.csi",
    output:
        bcf="called/{sample}_norm_qflt.bcf",
    params:
        qual_cutoff=get_qual_cutoff,
    # conda:
    #     "envs/bcftools.yaml"
    shell:
        # "vembrane tag -t low_qual='QUAL >= {params.qual_cutoff}' -O bcf {input.bcf}"
        "bcftools filter -i 'QUAL>{params.qual_cutoff}' -g 5 -Ob -o {output.bcf} {input.bcf}"


rule filter_het:
    input:
        bcf="called/{sample}_norm_qflt.bcf",
        csi="called/{sample}_norm_qflt.bcf.csi",
    output:
        bcf="called/{sample}_norm_qflt_het.bcf",
    # conda:
    #     "envs/bcftools.yaml"
    shell:
        "bcftools filter -i 'GT=\"hom\"' -Ob -o {output} {input.bcf}"


rule format_bcf:
    input:
        bcf=get_final_bcf,
        csi=lambda w: get_final_bcf(w, csi=True)
    output:
        "tsvs/{sample}.tsv",
    params:
        format=get_query_format,
    # conda:
    #     "envs/bcftools.yaml"
    shell:
        "bcftools query -f '{params.format}' -o {output} {input.bcf}"


rule make_consensus_genome:
    input:
        bcf=get_final_bcf,
        bcf_idx=lambda w: get_final_bcf(w, csi=True),
        ref=get_ref,
        ref_idx=lambda w: get_ref(w, fai=True),
    output:
        "seqs/{sample}.fa",
    # conda:
    #     "envs/bcftools.yaml"
    # params:
        # flt="FILTER!=\".\""
    shell:
        "bcftools consensus -f {input.ref} {input.bcf} -o {output}"


rule merge_tsvs:
    input:
        get_tsvs_to_merge,
        expand("tsvs/{sample}.tsv", sample=samples),
    output:
        "tsvs/merged.tsv",
    shell:
        "cat {input} > {output}"
