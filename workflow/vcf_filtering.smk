rule select_biallelic:
    input:
        "".join(["data/calls/{group}_", f"{config['calling']['caller']}_raw.bcf"]),
    output:
        "data/calls/{group}_{caller}_{suffix}_biallelic.bcf",
    threads: 4
    shell:
        "bcftools view -m2 -M2 -Ob -o {output} {input}"


rule select_snps:
    input:
        get_call_file,
    output:
        "data/calls/{group}_{caller}_{suffix}_snps.bcf",
        # "".join(["{input.calls}".strip(".bcf"), "_snps.bcf"]),
    threads: 4
    shell:
        "bcftools view -v snps -Ob -o {output} {input}"


# rule bcftools_norm:
#     input:
#         calls=get_call_file,
#         ref=get_ref,
#     output:
#         "data/calls/{group}_{caller}_{suffix}_norm.bcf",
#     resources:
#         time="2:00:00",
#     shell:
#         "bcftools norm -f {input.ref} -m '-' -Ob -o {output} {input.calls}"


rule filter_low_quality:
    input:
        calls=get_call_file,
        csi=lambda w: get_call_file(w, csi=True),
    output:
        "data/calls/{group}_{caller}_{suffix}_flt.bcf",
        # "".join(["{input.calls}".strip(".bcf"), "_flt.bcf"]),
    params:
        qual_cutoff=lambda w: get_filtering_criteria(w),
    shell:
        "bcftools filter -i {params.qual_cutoff} -s LowQual -g 5 -Ob -o {output} {input.calls}"


rule filter_homozygous:
    input:
        calls=get_call_file,
        csi=lambda w: get_call_file(w, csi=True),
    output:
        "data/calls/{group}_{caller}_{suffix}_homozygous.bcf",
        # "".join(["{input.calls}".strip(".bcf"), "_homozygous.bcf"]),
    shell:
        "bcftools filter -i 'GT=\"1/1\"' -Ob -o {output} {input.calls}"


# rule format_bcf:
#     input:
#         # bcf=get_final_bcf,
#         # csi=lambda w: get_final_bcf(w, csi=True),
#         bcf="data/calls/{comparison}/{condition}.vcf.gz",
#         csi="data/calls/{comparison}/{condition}.vcf.gz.csi",
#     output:
#         # "data/tsvs/{sample}.tsv",
#         "data/tsvs/{comparison}_{condition}.tsv",
#     params:
#         format=get_query_format,
#     shell:
#         "bcftools query -f '{params.format}' -o {output} {input.bcf}"

rule format_bcf:
    input:
        bcf=get_bcf_wildcard,
        csi=lambda w: get_bcf_wildcard(w, csi=True),
        # bcf="data/calls/{comparison}/{condition}.vcf.gz",
        # csi="data/calls/{comparison}/{condition}.vcf.gz.csi",
    output:
        "data/tsvs/{group}.tsv",
        # "data/tsvs/{comparison}_{condition}.tsv",
    params:
        format=get_query_format,
    shell:
        "bcftools query -f '{params.format}' -o {output} {input.bcf}"

# rule separate_into_samples:
#     input:
#         get_group_from_sample,
#     output:
#         "data/calls/{sample}_{suffix}.bcf",
#     resources:
#         time="2:00:00",
#     shell:
#         "bcftools view -s {wildcards.sample} -a -o {output.bcfs} {input}"


# rule compare_snps:
#     input:
#         bcf_list=get_comparisons,
#         bcf_idxs=lambda w: get_comparisons(w, csi=True),
#     output:
#         expand(
#             "data/calls/comparisons/{{comparison}}/{{f1}}_{condition}_{{f2}}.bcf",
#             condition=["not_in", "in"],
#             ),
#         expand(
#             "data/calls/comparisons/{{comparison}}/{{f2}}_{condition}_{{f1}}.bcf",
#             condition=["not_in", "in"],
#             )
#     params:
#         basenames=lambda w: get_comparisons(w, basename=True),
#     shell:
#         "scripts/compare_snps.sh {input.bcf_list} {wildcards.comparison} {params.basenames}"

# rule append_suffixes_to_comparisons:
#     input:
#         "data/calls/comparisons/{comparison}/{f1}_{condition}_{f2}.bcf",
#     output:
#         "data/calls/comparisons/{comparison}/{f1}_{condition}_{f2}_{suffix}.bcf",
#     params:
#         sfix=lambda w: get_suffix_for_comparison(w)
#     shell:
#         "mv {input} {input}_{params.sfix}.bcf"

# ruleorder: bcftools_call > separate_into_samples
# ruleorder: separate_into_samples > append_suffixes_to_comparisons
# ruleorder: append_suffixes_to_comparisons > select_biallelic


rule make_consensus_genome:
    input:
        bcf=lambda w: get_final_bcf(w.sample),
        bcf_idx=lambda w: get_final_bcf(w.sample, csi=True),
        ref=get_ref,
        ref_idx=lambda w: get_ref(w, fai=True),
    output:
        "data/consensus/{sample}.fa",
    shell:
        "bcftools consensus -f {input.ref} {input.bcf} -o {output}"


rule merge_tsvs:
    input:
        get_tsvs_to_merge,
    output:
        "data/tsvs/merged.tsv",
    shell:
        "cat {input} > {output}"
