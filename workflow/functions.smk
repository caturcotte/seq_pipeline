import pandas as pd
import numpy as np
import re
from pathlib import Path

sample_sheet = pd.read_csv("sample_sheet.csv")

# get sample field information to create sample names
sample_sheet['sample_num'] = round(sample_sheet['sample_num'])
sample_fields = list(zip(
    sample_sheet["condition"],
    sample_sheet["sample_type"].fillna(''),
    sample_sheet["sample_num"].fillna(''),
))
s_fields = []
for i in sample_fields:
    s_fields.append([j for j in i if j])
for i in s_fields:
    if len(i) > 2:
        i[2] = str(int(i[2])).zfill(3)
sample_names = ['-'.join(i) for i in s_fields]
sample_sheet['sample'] = sample_names

# identify ref and alt parent, make sure they exist as samples
ref_parent_name = config["ref_parent"]
alt_parent_name = config["alt_parent"]
parent_names = [ref_parent_name, alt_parent_name]
try:
    ref_parent = sample_sheet.loc[sample_sheet["sample"] == ref_parent_name]
    alt_parent = sample_sheet.loc[sample_sheet["sample"] == alt_parent_name]
except:
    raise OSError("One or more parents were not found in the sample sheet.")

samples = list(sample_sheet["sample"])
progeny = [i for i in samples if i not in parent_names]


# list all lane IDs for one sample
def list_all_idens(sample_sheet):
    prefixes = []
    for sample_name in list(sample_sheet["sample"]):
        [prefixes.append(i) for i in get_ids_for_sample(sample_name)]
    return prefixes


# return sample row given sample name
def get_sample(sample_name):
    return sample_sheet.loc[sample_sheet["sample"] == sample_name]


# return split fields from sample name
def get_sample_fields(sample_name):
    fields = sample_name.split("-")
    return dict(zip(["condition", "sample_type", "sample_num"], fields))


# check if reads are illumina or nanopore
def is_short_read(sample_name):
    sample = get_sample(sample_name)
    platform = sample["platform"].iloc[0].lower()
    if platform == "illumina":
        return True
    elif platform == "nanopore":
        return False
    else:
        raise ValueError(
            f"Invalid platform specified for sample {sample} in sample_sheet.csv. Platform must be either Illumina or Nanopore."
        )


# check if a sample has lane IDs or not
def has_lane_id(sample_name):
    sample = get_sample(sample_name)
    files = find_all_fastqs(sample_name)
    ids = get_ids_for_sample(sample_name)
    return bool(ids)


# get reference file name/base name, non input function, name depends on whether repeats are masked
def get_ref_no_input(base=False, masked=config['mask_repeats'], fai=False):
    name = ["data/resources/", config["ref_name"]]
    if masked:
        name.append("_masked")
    if base:
        return "".join(name)
    name.append(".fa")
    if not fai:
        return "".join(name)
    name.append(".fai")
    return "".join(name)


# get reference file name/base name, used as input function so needs to have wildcards
def get_ref(w, base=False, masked=config["mask_repeats"], fai=False):
    return get_ref_no_input(base, masked, fai)


# return the path at which a sample's reads are located
def get_data_path(sample_name):
    sample = get_sample(sample_name)
    location = sample["location"].iloc[0]
    if location not in config["data_locations"]:
        raise ValueError(f"Sample location {location} not found in config file.")
    filepath = Path(config["data_locations"][location])
    if not filepath.exists():
        raise OSError(f"File directory {str(filepath)} not found.")
    dirs = [f"usftp21.novogene.com/01.RawData/{sample_name}", sample_name]
    barcode = sample["barcode"]
    if barcode.notnull().all():
        barcode = int(barcode.iloc[0])
        dirs.append(f"fastq_pass/{barcode:02}")
    # else:
    #     path_for_sample = 
    paths = [filepath.joinpath(i) for i in dirs]
    valid_file_paths = [i for i in paths if i.exists()]
    if len(valid_file_paths) > 1:
        raise OSError(f"Multiple conflicting filepaths found in {location}.")
    elif valid_file_paths == []:
        raise OSError(f"Directories containing reads not found in {location}.")
    return valid_file_paths[0]


# search each file in a directory for a regex match
def regex_over_dir(path, regex):
    p = re.compile(regex)
    matches = []
    for i in path.iterdir():
        match = p.search(str(i))
        if match is not None:
            matches.append(match)
    return matches


# find all fastqs associated with a sample
def find_all_fastqs(sample_name):
    path = get_data_path(sample_name)
    sample = get_sample(sample_name)
    # if sample["barcode"].notnull().all():
    #     path.rglob(barcode
    #     regex = rf"(.*{sample_name}_\w*(?:_[12])?\.f(?:ast)?q(?:\.gz)?)"
    # else:Vk
    regex = rf"(.*{sample_name}_[a-zA-Z0-9_-]*(?:_[12])?\.f(?:ast)?q(?:\.gz)?)"
    matches = regex_over_dir(path, regex)
    if not matches:
        raise OSError(f"No reads found in {str(path)} for sample {sample_name}.")
    return [i[0] for i in matches]


def get_ids_for_sample(sample_name):
    path = get_data_path(sample_name)
    regex = rf"(?<={sample_name}_)([a-zA-Z0-9_-]*)(?=_[12]\.f(?:ast)?q(?:\.gz)?)"
    matches = regex_over_dir(path, regex)
    if matches:
        return [i[0] for i in matches]
    else:
        return [sample_name]


def find_fastq(path):
    options = [".fq", ".fq.gz", ".fastq", ".fastq.gz"]
    for i in options:
        file = path.with_suffix(i)
        if file.is_file():
            return str(file)
    raise OSError(f"No fastq file found for {path}.")


# get the format for transforming vcfs into tsvs
def get_query_format(w):
    cols = ["%SAMPLE", "%CHROM", "%POS", "%REF", "%ALT", "%QUAL", "%GT", "%DP", "%AD"]
    query = "\t".join(cols)
    return "[" + query + "]\n"


"""
Input functions
"""


def find_all_fastqs_for_id(w):
    all_fastqs = find_all_fastqs(w.sample)
    if w.sample == w.iden:
        return all_fastqs
    else:
        regex = rf"(.*{w.sample}_{w.iden}{w.read}\.f(?:ast)?q(?:\.gz)?)"
        path = get_data_path(w.sample)
        matches = regex_over_dir(path, regex)
        return [i[0] for i in matches]


def get_alns_to_merge(w):
    return [
        str(Path("data/alignments/").joinpath(f"{w.sample}_{i}_sort.bam"))
        for i in get_ids_for_sample(w.sample)
    ]


def get_fastqc_files(w):
    if is_short_read(w.sample):
        return expand(
            "data/qc/fastqc/{sample}_{read}_fastqc.zip", sample=samples, read=[1, 2]
        )
    else:
        return expand("data/qc/fastqc/{sample}_fastqc.zip", sample=samples)


# get required names of aligned reads, determines which aligner is used
def get_aligned_reads(w):
    if is_short_read(w.sample):
        return f"data/alignments/{w.sample}_{w.iden}_bt2.bam"
    else:
        return f"data/alignments/{w.sample}_{w.iden}_mm2.bam"


def get_final_bams(w, bai=False):
    sample = get_sample(w.sample)
    file = ["data/alignments/{sample}"]
    if sample["PCR"]:
        file.append("_markdup")
    file.append(".bam")
    if bai:
        file.append(".bai")
    return "".join(file)


def get_reads_to_map(w):
    sample = get_sample(w.sample)
    if has_lane_id(w.sample):
        f = "data/reads/{sample}_{iden}"
    else:
        f = "data/reads/{sample}_{sample}"
    if is_short_read(w.sample):
        basename = [f"{f}_1", f"{f}_2"]
    else:
        basename = [f]
    final = ["".join([i, ".fq.gz"]) for i in basename]
    return final


# get ref index for bowtie2
def get_ref_bowtie2(w):
    return multiext(
        get_ref(w, base=True), ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"
    )


# get ref index for minimap2
def get_ref_minimap2(w):
    return "".join([get_ref(w), ".mmi"])
