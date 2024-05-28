import os
from itertools import combinations, chain
from snakemake.utils import validate


configfile: "config.yaml"
validate(config, "schemas/config.schema.yaml")


include: "workflow/functions.smk"
include: "workflow/alignment.smk"
include: "workflow/calling.smk"
include: "workflow/misc.smk"
include: "workflow/qc.smk"
include: "workflow/vcf_filtering.smk"

validate(sample_sheet, "schemas/sample.schema.yaml")

localrules:
    all,
    mv_masked_ref,
    mv_nolane_bams,
    symlink_fqs_single,
    symlink_ref,
    symlink_fqs_paired,
    mv_fastas,


prefixes = []
for sample_name in list(sample_sheet["sample"]):
    [prefixes.append(i) for i in get_ids_for_sample(sample_name)]

samples_and_ref = samples + [config["ref_name"]]

# samples_or_groups = []
# for i in config["analysis"]["comparisons"]:
#     for j in config["analysis"]["comparisons"][i]["datasets"]:
#         samples_or_groups.append(j)
# samples_and_groups = list(set(list(groups.keys()) + samples_or_groups))

wildcard_constraints:
    group="|".join([i for i in groups.keys()]),
    sample="|".join([i for i in samples]),
    ref=get_ref(base=True),
    iden="|".join(prefixes),
    sample_or_ref="|".join(samples_and_ref),
    suffix="|".join(suffix_constraints),
    caller=config['calling']['caller'],
    # sample_or_group="|".join(samples_and_groups),
    i=r"\d+",
    label="|".join(get_labels()),
    # chrom="|".join([i for i in chroms]),
    # ref="|".join([i for i in ref_names]),


rule all:
    input:
        "data/qc/multiqc.html",
        get_final_output,
