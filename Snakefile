import os


configfile: "config.yaml"


include: "workflow/functions.smk"
include: "workflow/alignment.smk"
include: "workflow/calling.smk"
include: "workflow/misc.smk"
include: "workflow/vcf_filtering.smk"


localrules:
    all,
    mv_masked_ref,
    symlink_fqs_single,
    symlink_ref,
    symlink_fqs_paired,


wildcard_constraints:
    group="|".join([i for i in groups.keys()]),
    sample="|".join([i for i in samples]),
    ref=config['ref_name']
    # chrom="|".join([i for i in chroms]),
    # ref="|".join([i for i in ref_names]),


rule all:
    input:
        get_final_output,
