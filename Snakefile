import os


configfile: "config.yaml"


include: "workflow/functions.smk"
include: "workflow/alignment.smk"


if config["caller"] == "freebayes" or config["bcftools_opts"]["call_as_groups"]:

    include: "workflow/calling_group.smk"

else:

    include: "workflow/calling_single.smk"


include: "workflow/misc.smk"
include: "workflow/vcf_filtering.smk"


if config["ndj_analysis"]:

    include: "workflow/ndj_analysis.smk"

localrules:
    all,
    symlink_fqs_single,
    symlink_ref,
    symlink_fqs_paired

wildcard_constraints:
    group="|".join([i for i in groups.keys()]),
    sample="|".join([i for i in samples]),
    chrom="|".join([i for i in chroms]),
    ref="|".join([i for i in ref_names]),
    masked="|".join(["masked", ""]),


rule all:
    input:
        get_final_output,
