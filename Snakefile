import os


configfile: "config.yaml"


include: "workflow/functions.smk"
include: "workflow/alignment.smk"
include: "workflow/calling.smk"
include: "workflow/misc.smk"
include: "workflow/qc.smk"
include: "workflow/vcf_filtering.smk"


localrules:
    all,
    mv_masked_ref,
    symlink_fqs_single,
    symlink_ref,
    symlink_fqs_paired,
    mv_fastas,


prefixes = []
for sample_name in list(sample_sheet["sample"]):
    [prefixes.append(i) for i in get_ids_for_sample(sample_name)]

samples_and_ref = samples + [config["ref_name"]]


wildcard_constraints:
    group="|".join([i for i in groups.keys()]),
    sample="|".join([i for i in samples]),
    ref=get_ref_basename(),
    iden="|".join(prefixes),
    sample_or_ref="|".join(samples_and_ref),
    i=r"\d+",
    label="|".join(get_labels()),
    # chrom="|".join([i for i in chroms]),
    # ref="|".join([i for i in ref_names]),


rule all:
    input:
        "data/qc/multiqc.html",
        get_final_output,
