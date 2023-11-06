ref_parent = config["ndj_analysis_opts"]["parents"]["ref_parent"]


rule bcftools_call_for_ref_parent:
    input:
        alns=f"mapped/{ref_parent}_sort_dedup.bam",
        aln_idxs=f"mapped/{ref_parent}_sort_dedup.bam.bai",
        ref=f"resources/{config['ref_name']}.fa",
        ref_idx=f"resources/{config['ref_name']}.fa.fai",
    output:
        f"called/ref_parent/{ref_parent}_unprocessed.bcf",
    conda:
        "envs/bcftools.yaml"
    threads: 16
    shell:
        "bcftools call -cv -f {input.ref} {input.alns} -o {output} --threads {threads}"


rule filter_ref_parent:
    input:
        bcf=f"called/ref_parent/{ref_parent}_unprocessed.bcf",
        ref=f"resources/{config['ref_name']}.fa",
        ref_idx=f"resources/{config['ref_name']}.fa.fai",
    output:
        bcf=f"called/ref_parent/{ref_parent}_norm_flt.bcf",
        csi=f"called/ref_parent/{ref_parent}_norm_flt.bcf.csi",
    resources:
        time="2:00:00",
    params:
        extra="--write-index",
        qual_cutoff=lambda w: get_qual_cutoff(w),
    conda:
        "envs/bcftools.yaml"
    shell:
        "bcftools norm -ad all -f {input.ref} {input.bcf} | bcftools filter -i 'QUAL < {params.qual_cutoff}' -s lowQual -o {output.bcf} --write-index -"


rule make_ref_parent_genome:
    input:
        bcf=f"called/ref_parent/{ref_parent}_norm_flt.bcf",
        csi=f"called/ref_parent/{ref_parent}_norm_flt.bcf.csi",
        ref=f"resources/{config['ref_name']}.fa",
        ref_idx=f"resources/{config['ref_name']}.fa.fai",
    output:
        f"resources/{ref_parent}.fa",
    conda:
        "envs/bcftools.yaml"
    shell:
        "bcftools consensus -f {input.ref} {input.bcf} -e 'FILTER != .' -o {output}"


rule get_unique_snps_progeny:
    input:
        progeny=get_final_bcf,
        idx=lambda w: get_final_bcf(w, csi=True),
        alt_parent="called/"
        + config["ndj_analysis_opts"]["parents"]["alt_parent"]
        + "_norm_qflt_het.bcf",
    output:
        expand("called/isecs/{{sample}}/000{num}.vcf", num=[0, 1, 2, 3]),
    conda:
        "envs/bcftools.yaml"
    shell:
        "bcftools isec {input.progeny} {input.alt_parent} -p called/isecs/{wildcards.sample}"


rule format_progeny_common_snp_vcfs:
    input:
        "called/isecs/{sample}/0000.vcf",
    output:
        "tsvs/{sample}_common.tsv",
    params:
        format=get_query_format,
    conda:
        "envs/bcftools.yaml"
    shell:
        "bcftools query -f {params.format} -o {output}"


rule merge_progeny_common_snp_vcfs:
    input:
        get_progeny_tsvs,
    output:
        "tsvs/progeny_common.tsv",
    conda:
        "envs/bcftools.yaml"
    params:
        format=get_query_format,
    shell:
        "bcftools query -f {params.format} -o {output}"
