import pandas as pd
import numpy as np
import re
from pathlib import Path


# get reference file name/base name, used as input function so needs to have wildcards
def get_ref(w, config=config, base=False, fai=False):
    return get_ref_no_input(config, base, fai)


def get_all_pod5_files(w, config=config):
    in_dir = Path(config['base_dir'])
    pod5s = in_dir.glob('*/*/pod5/*.pod5')
    return [str(i) for i in pod5s]

def get_all_pod5_tags(w, config=config):
    in_dir = Path(config['base_dir'])
    pod5s = in_dir.glob('*/*/pod5/*.pod5')
    return [str(i.stem) for i in pod5s]

def find_pod5(w, config=config):
    pod5s=get_all_pod5s(config)
    pod5 = [i for i in pod5s if w.pod5 in i]
    if len(pod5) > 1:
        raise ValueError('Duplicate pod5s found.')
    return pod5[0]

def get_progeny_run_folder_combos(w, sample_sheet=sample_sheet):
    files = []
    for i in list(sample_sheet['sample']):
        run_folder = sample_sheet.loc[sample_sheet['sample'] == i, 'ont_folder'].iloc[0]
        files.append(f"data/parquets/{run_folder}/{i}_rearranged.parquet")
    return files
        # expand("data/parquets/{run_folder}/{progeny}_rearranged.parquet", run_folder=folders, progeny=all_progeny, allow_missing=True),

# def get_folder_from_pod5(prefix, config, df, folders):
#     pod5_dict = get_pod5_dict(folders, config, df)
#     for i in pod5_dict:
#         if w.prefix in pod5_dict[i].stem:
#             return i
#     raise ValueError(f'{prefix} not found in any folder.')

# def get_pod5_file(w, config=config, sample_sheet=sample_sheet):
#     pod5 = Path(config['ont_base_location']).joinpath(
#         config['experiment_id'],
#         get_sample_id_from_folder(sample_sheet, w.folder),
#         w.folder,
#         f"{w.pod5_name}.pod5"
#         )
#     return str(pod5)


def get_pod5_with_tag(w):
    pod5s = get_all_pod5s(w.run_info)
    prefix = [i for i in pod5s if w.tag in i]
    if len(prefix) > 1:
        raise ValueError("Non-unique pod5s found")
    return prefix[0]
    # pod5_names = pod5_dict[w.folder]
    # return expand(
    #     'data/pod5s/{{folder}}/{pod5_name}.pod5',
    #     pod5_name=pod5_names
    #     )
    # sample_ids = get_all_sample_ids(sample_sheet)
    # sample_dirs = [exp_dir.joinpath(i) for i in sample_ids]
    # for sample_dir in sample_dirs:
    #     if not sample_dir.is_dir():
    #         raise OSError(f'Directory for sample {str(sample_dir)} not found in {str(exp_dir)}.')
    #     # get every folder that's in that sample directory
    #     folders = list(sample_sheet.loc[sample_sheet['sample_id'].isin(sample_dir.stem), 'ont_folder'].unique())


def get_pod5_files(w):
    (prefixes,) = glob_wildcards(f"data/pod5s/ndj_progeny/{w.folder}/{{prefix}}.pod5")
    return expand("data/pod5s/ndj_progeny/{{folder}}/{prefix}.pod5", prefix=prefixes)


# get the format for transforming vcfs into tsvs
def get_query_format(w, sample_sheet=sample_sheet, config=config):
    cols = ["%SAMPLE", "%CHROM", "%POS", "%REF", "%ALT", "%QUAL", "%GT", "%DP", "%AD"]
    if not is_short_read(w.sample, sample_sheet):
        cols.append("%PS")
    query = "\t".join(cols)
    return "[" + query + "]\n"


def find_all_fastqs_for_id(w, sample_sheet=sample_sheet, config=config):
    path = get_data_path(w.sample, sample_sheet, config)
    all_fastqs = get_fastqs_with_gz(w.sample, sample_sheet, config)
    return [path.joinpath(i) for i in all_fastqs]
    # if w.sample == w.iden:
    return all_fastqs
    # else:
    #     # regex = rf"(.*{w.sample}_{w.iden}{w.read}\.f(?:ast)?q(?:\.gz)?)"
    #     regex = rf"(.*{w.sample}_{w.iden}{w.read}\.f(?:ast)?q\.gz)"
    #     path = get_data_path(w.sample, sample_sheet, config)
    #     matches = regex_over_dir(path, regex)
    #     return [i[0] for i in matches]


def get_fastq_from_id(w, sample_sheet=sample_sheet, config=config):
    return get_fastq(w.sample, w.iden, w.read, sample_sheet, config)


def get_sedex(w, sample_sheet=sample_sheet):
    if is_short_read(w.sample, sample_sheet):
        return "'1s/^/sample\\tchromosome\\tposition\\treference\\tvariant\\tquality\\tgenotype\\tdepth\\tallele_depth\\n/'"
    else:
        return "'1s/^/sample\\tchromosome\\tposition\\treference\\tvariant\\tquality\\tgenotype\\tdepth\\tallele_depth\\tphase_set\\n/'"


def check_if_short_read(w, sample_sheet=sample_sheet):
    if w.sample in list(sample_sheet["sample"]):
        return False
    return True


#     if is_short_read(w.sample, sample_sheet):
#         return True
#     return False
# def get_alns_to_merge(w, sample_sheet=sample_sheet, config=config):
#     return [
#         str(Path("data/alignments/").joinpath(f"{w.sample}_{i}_sort.bam"))
#         for i in get_ids_for_sample(w.sample, sample_sheet, config)
# ]


def get_fastqc_files(w, sample_sheet=sample_sheet):
    if is_short_read(w.sample, sample_sheet):
        return expand(
            "data/qc/fastqc/{sample}_{read}_fastqc.zip", sample=samples, read=[1, 2]
        )
    else:
        return expand("data/qc/fastqc/{sample}_fastqc.zip", sample=samples)


# get required names of aligned reads, determines which aligner is used
def get_aligned_reads(w, sample_sheet=sample_sheet):
    if is_short_read(w.sample, sample_sheet):
        return f"data/alignments/{w.sample}_{w.iden}_bt2.bam"
    else:
        return f"data/alignments/{w.sample}_{w.iden}_mm2.bam"


def get_final_bams(w, bai=False, sample_sheet=sample_sheet):
    sample = get_sample(w.sample, sample_sheet)
    file = ["data/alignments/{sample}"]
    if sample["PCR"]:
        file.append("_markdup")
    file.append(".bam")
    if bai:
        file.append(".bai")
    return "".join(file)


# updated to get bam from dorado instead of fastq
def get_reads_to_map(w, sample_sheet=sample_sheet, config=config):
    sample = get_sample(w.sample, sample_sheet)
    if has_lane_id(w.sample, sample_sheet, config):
        f = "data/reads/{sample}_{iden}"
    else:
        f = "data/reads/{sample}_{sample}"
    if is_short_read(w.sample, sample_sheet):
        basename = [f"{f}_1", f"{f}_2"]
        # return ["".join([i, ".fq.gz"]) for i in basename]
    else:
        # return "".join([f, ".bam"])
        basename = [f]
    final = ["".join([i, ".fq.gz"]) for i in basename]
    return final


# get ref index for bowtie2
def get_ref_bowtie2(w):
    return multiext(
        get_ref(w, base=True),
        ".1.bt2",
        ".2.bt2",
        ".3.bt2",
        ".4.bt2",
        ".rev.1.bt2",
        ".rev.2.bt2",
    )


# get ref index for minimap2
def get_ref_minimap2(w):
    return "".join([get_ref(w), ".mmi"])


def get_bcf_file(w, sample_sheet=sample_sheet, csi=False):
    if is_short_read(w.sample, sample_sheet):
        caller = "bcftools"
        return f"data/calls/{w.sample}_raw_{caller}.bcf"
    else:
        caller = "phased"
    if csi:
        return f"data/calls/{w.sample}_raw_{caller}.bcf.csi"
    return f"data/calls/{w.sample}_raw_{caller}.bcf"


def aggregate_input(w, breaks=False):
    # decision based on content of output file
    # Important: use the method open() of the returned file!
    # This way, Snakemake is able to automatically download the file if it is generated in
    # a cloud environment without a shared filesystem.
    with open(checkpoints.make_mosdepth_df.get().output[0], "r") as file:
        smps = [i.rstrip() for i in file]
    passed_samples = []
    for i in smps:
        if breaks:
            passed_samples.append(
                f"data/tiger/{i}/{i}.CO_estimates.smoothed.breaks.txt"
            )
        else:
            passed_samples.append(f"data/tiger/{i}/{i}.CO_estimates.txt")
    return passed_samples
