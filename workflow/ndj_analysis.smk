ref_parent = config["ndj_analysis_opts"]["parents"]["ref_parent"]
rule bcftools_mpileup_ref_parent:
    input:
        alns=f"data/alignments/{ref_parent}_sort_dedup.bam",
        aln_idxs=f"data/alignments/{ref_parent}_sort_dedup.bam.bai",
        ref=get_ref_for_ref_parent,
        ref_idx=lambda w: get_ref_for_ref_parent(w, fai=True),
    output:
        f"data/calls/ref_parent/{ref_parent}_bcftools_pileup.bcf",
    params:
        regions=lambda w: get_regions_to_call(w),
        extra="--max-depth 200 --min-BQ 15"
    threads: 16
    # conda:
    #     "envs/bcftools.yaml"
    shell:
        "bcftools mpileup -f {input.ref} -r {params.regions} {params.extra} -Ou -o {output} --threads {threads} {input.alns}"


rule bcftools_call_for_ref_parent:
    input:
        bcf=f"data/calls/{ref_parent}_bcftools_pileup.bcf",
        txt=f'.tmp/single_call_{ref_parent}.txt'
    output:
        f"data/calls/ref_parent/{ref_parent}_unprocessed.bcf",
    # conda:
    #     "envs/bcftools.yaml"
    params:
        regions=lambda w: get_regions_to_call(w),
        extra="--max-depth 200 --min-BQ 15"
    threads: 16
    shell:
        "bcftools call -cv -o {output} --threads {threads} {input.bcf}"

rule bcftools_norm_ref_parent:
    input:
        calls=f"data/calls/ref_parent/{ref_parent}_unprocessed.bcf",
        ref=get_ref_for_ref_parent,
    output:
        bcf=temp(f"data/calls/ref_parent/{ref_parent}_norm.bcf"),
    # conda:
    #     "envs/bcftools.yaml"
    resources:
        time="2:00:00",
    shell:
        "bcftools norm -f {input.ref} -m '-' -Ob -o {output.bcf} {input.calls}"

rule filter_ref_parent:
    input:
        bcf=f"data/calls/ref_parent/{ref_parent}_norm.bcf",
        idx=f"data/calls/ref_parent/{ref_parent}_norm.bcf.csi",
        # ref=f"resources/{config['ref_name']}.fa",
        # ref_idx=f"resources/{config['ref_name']}.fa.fai",
    output:
        bcf=f"data/calls/ref_parent/{ref_parent}_norm_flt.bcf",
    resources:
        time="2:00:00",
    params:
        qual_cutoff=lambda w: get_qual_cutoff(w),
    # conda:
    #     "envs/bcftools.yaml"
    shell:
        "bcftools filter -i 'QUAL>{params.qual_cutoff} & GT=\"hom\"' -g 5 -Ob -o {output} {input.bcf}"
        # "vembrane tag -t low_qual='QUAL >= {params.qual_cutoff}' -t het='is_het' -o {output} -O bcf {input.bcf}"
        # "bcftools filter -e'QUAL > {params.qual_cutoff} & GT=\"1/1\"' -g 5 -s lowQual -Ob -o {output.bcf} {input.bcf}"


rule make_ref_parent_genome:
    input:
        bcf=f"data/calls/ref_parent/{ref_parent}_norm_flt.bcf",
        csi=f"data/calls/ref_parent/{ref_parent}_norm_flt.bcf.csi",
        ref=get_ref_for_ref_parent,
        ref_idx=lambda w: get_ref_for_ref_parent(w, fai=True),
    output:
        f"resources/{ref_parent}.fa",
    # conda:
    #     "envs/bcftools.yaml"
    shell:
        "bcftools consensus -f {input.ref} {input.bcf} -o {output}"


rule get_unique_snps_progeny:
    input:
        progeny=get_final_bcf,
        idx=lambda w: get_final_bcf(w, csi=True),
        alt_parent="data/calls/"
        + config["ndj_analysis_opts"]["parents"]["alt_parent"]
        + "_norm_qflt_het.bcf",
        alt_parent_idx="data/calls/"
        + config["ndj_analysis_opts"]["parents"]["alt_parent"]
        + "_norm_qflt_het.bcf.csi",
    output:
        expand("data/calls/isecs/{{sample}}/000{num}.vcf", num=[0, 1, 2, 3]),
    # conda:
    #     "envs/bcftools.yaml"
    shell:
        "bcftools isec {input.progeny} {input.alt_parent} -p data/calls/isecs/{wildcards.sample}"


rule format_progeny_common_snp_vcfs:
    input:
        "data/calls/isecs/{sample}/0000.vcf",
    output:
        "data/tsvs/{sample}_common.tsv",
    params:
        format=get_query_format,
    # conda:
    #     "envs/bcftools.yaml"
    shell:
        "bcftools query -f '{params.format}' -o {output} {input}"


rule merge_progeny_common_snp_vcfs:
    input:
        get_progeny_tsvs,
    output:
        "data/tsvs/progeny_common.tsv",
    # conda:
    #     "envs/bcftools.yaml"
    params:
        format=get_query_format,
    shell:
        "cat {input} > {output}"
