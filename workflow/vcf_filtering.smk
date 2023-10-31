rule bcftools_norm:
    input:
        calls=get_caller,
        ref=get_ref,
    output:
        bcf=temp("called/{sample}_norm.bcf"),
        csi=temp("called/{sample}_norm.bcf.csi"),
    resources:
        time="2:00:00",
    params:
        extra="--write-index",
    wrapper:
        "v2.3.1/bio/bcftools/norm"


rule filter_low_quality:
    input:
        bcf="called/{sample}_norm.bcf",
        csi="called/{sample}_norm.bcf.csi",
    output:
        bcf="called/{sample}_norm_qflt.bcf",
        csi="called/{sample}_norm_qflt.bcf.csi",
    params:
        qual_cutoff=get_qual_cutoff,
    shell:
        "bcftools filter -i 'QUAL < 150' -s lowQual -o {output.bcf} --write-index {input.bcf}"


rule filter_het:
    input:
        bcf="called/{sample}_norm_qflt.bcf",
        csi="called/{sample}_norm_qflt.bcf.csi",
    output:
        bcf="called/{sample}_norm_qflt_het.bcf",
        csi="called/{sample}_norm_qflt_het.bcf.csi",
    shell:
        "bcftools filter -i 'GT != 1/1' -s het -o {output.bcf} --write-index {input.bcf}"


rule format_bcf:
    input:
        get_final_bcf,
    output:
        "tsvs/{sample}.tsv",
    params:
        format=get_query_format,
    shell:
        "bcftools query -f {params.format} -o {output}"


rule make_consensus_genome:
    input:
        bcf=get_final_bcf,
        bcf_idx=lambda w: get_final_bcf(w, csi=True),
        ref=get_ref,
        ref_idx=lambda w: get_ref(w, fai=True),
    output:
        "seqs/{sample}.fa",
    shell:
        "bcftools consensus -f {input.ref} {input.bcf} -e 'FILTER != .' -o {output}"


rule merge_tsvs:
    input:
        get_tsvs_to_merge,
        expand("tsvs/{sample}.tsv", sample=samples),
    output:
        "tsvs/merged.tsv",
    shell:
        "cat {input} > {output}"
