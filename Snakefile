import os


configfile: "config.yaml"


include: "workflow/functions.smk"
include: "workflow/alignment.smk"
if config['call_as_groups']:
    include: "workflow/calling_group.smk"
else:
    include: "workflow/calling_single.smk"
include: "workflow/misc.smk"
include: "workflow/ndj_analysis.smk"
include: "workflow/vcf_filtering.smk"


wildcard_constraints:
    group="|".join(["parents", "progeny"]),
    sample="|".join([i for i in samples]),
    # parent="|".join([i for i in parents]),
    # progeny="|".join([i for i in progeny]),
    chrom="|".join([i for i in chroms]),
    ref="|".join([i for i in ref_names]),
    # name="|".join([j for i in groups for j in groups[i]]),


rule all:
    input:
        get_final_output,

