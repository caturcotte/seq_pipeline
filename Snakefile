import os
from itertools import combinations, chain
from snakemake.utils import validate


configfile: "config.yaml"


validate(config, "schemas/config.schema.yaml")


include: "workflow/functions.smk"
include: "workflow/alignment.smk"
include: "workflow/calling.smk"
include: "workflow/misc.smk"
include: "workflow/vcf_filtering.smk"


validate(sample_sheet, "schemas/sample.schema.yaml")


localrules:
    all,
    mv_masked_ref,
    mv_nolane_bams,
    add_header_to_tsv,
    symlink_ref,


wildcard_constraints:
    sample="|".join(list(sample_sheet["sample"])),
    aligner="|".join(["bt2", "mm2"]),
    ref=get_ref_no_input(base=True),
    iden="|".join(list_all_idens(sample_sheet)),
    i=r"\d+",


rule all:
    input:
        "data/tiger/tiger_output.parquet",
