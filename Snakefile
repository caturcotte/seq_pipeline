import os
import pandas as pd
from itertools import combinations, chain
from snakemake.utils import validate


configfile: "config.yaml"
# validate(config, "schemas/config.schema.yaml")


include: "workflow/sample_sheet_functions.smk"


sample_sheet = pd.read_csv('sample_sheet.csv')
# validate (sample_sheet, 'schemas/sample.schema.yaml')
sample_sheet = concat_sample_names(sample_sheet)
sample_sheet.to_csv('test.csv')

all_progeny = get_progeny_names(sample_sheet)
parents = get_parent_names(config)
samples = all_progeny + list(parents.values())
folders = list(sample_sheet['ont_folder'].unique())
run_folders = [*folders, "parents"]

# pod5s = get_pod5_dict(folders, config, sample_sheet)
# progeny_dict = get_progeny_dict(sample_sheet)
os = detect_platform()[0]
arch = detect_platform()[1]

include: "workflow/input_functions.smk"
include: "workflow/misc.smk"
include: "workflow/alignment.smk"
include: "workflow/calling.smk"
include: "workflow/vcf_filtering.smk"
include: "workflow/plotting.smk"
# include: "workflow/tiger.smk"
# include: "workflow/qc.smk"


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
    progeny="|".join(all_progeny),
    parent="|".join(parents),
    run_folder="|".join(run_folders),
    #aligner="|".join(["bt2", "mm2"]),
    #caller="|".join(['clair3', 'bcftools', 'phased']),

# ont_samples = ["".join(["WT-F1-0", str(i)]) for i in list(range(17, 80))]
rule all:
    input:
        expand("data/plots/{progeny}_phased.csv", progeny=all_progeny)
        # expand("data/tiger/{sample}/{sample}_plot1.csv", sample=ont_samples),
        # expand("data/tiger/{sample}/{sample}_plot2.csv", sample=ont_samples),
        # expand("data/tiger/{sample}/{sample}_intvls.csv", sample=ont_samples),
        # expand("data/calls/{sample}_raw_phased.bcf", sample=ont_samples),
        # "data/qc/multiqc.html",
        # "data/tiger_output/hmm_states.csv",
        # "data/tiger_output/hmm_intervals.csv",
        # "data/tiger_output/chr_2_states.csv"
