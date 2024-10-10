import os
import pandas as pd
from itertools import combinations, chain
from snakemake.utils import validate


configfile: "config.yaml"
validate(config, "schemas/config.schema.yaml")


include: "workflow/sample_sheet_functions.smk"


sample_sheet = pd.read_csv('sample_sheet.csv')
validate (sample_sheet, 'schemas/sample.schema.yaml')
sample_sheet = concat_sample_names(sample_sheet)

samples = get_sample_names(sample_sheet)
progeny = get_progeny_names(sample_sheet, config)
parents = get_parent_names(sample_sheet, config)


include: "workflow/input_functions.smk"
include: "workflow/alignment.smk"
include: "workflow/calling.smk"
include: "workflow/misc.smk"
include: "workflow/vcf_filtering.smk"
include: "workflow/tiger.smk"
include: "workflow/qc.smk"


localrules:
    all,
    mv_masked_ref,
    mv_nolane_bams,
    add_header_to_tsv,
    symlink_ref,
    make_chrom_lengths_file,
    fix_break_files


wildcard_constraints:
    sample="|".join(samples),
    progeny="|".join(progeny),
    aligner="|".join(["bt2", "mm2"]),
    caller="|".join(['clair3', 'bcftools', 'phased']),
    ref=get_ref_no_input(config, base=True, fai=False),
    read="|".join(["_1", "_2", ""]),
    iden="|".join(list_all_idens(sample_sheet, config)),
    i=r"\d+",

ont_samples = ["".join(["WT-F1-0", str(i)]) for i in list(range(17, 80))]
rule all:
    input:
        expand("data/tiger/{sample}/{sample}_phased.csv", sample=progeny)
        # expand("data/tiger/{sample}/{sample}_plot1.csv", sample=ont_samples),
        # expand("data/tiger/{sample}/{sample}_plot2.csv", sample=ont_samples),
        # expand("data/tiger/{sample}/{sample}_intvls.csv", sample=ont_samples),
        # expand("data/calls/{sample}_raw_phased.bcf", sample=ont_samples),
        # "data/qc/multiqc.html",
        # "data/tiger_output/hmm_states.csv",
        # "data/tiger_output/hmm_intervals.csv",
        # "data/tiger_output/chr_2_states.csv"
