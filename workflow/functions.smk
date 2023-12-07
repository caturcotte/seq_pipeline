import pandas as pd
import numpy as np
import os
import re
import sys

print(sys.version)


reference = config["reference"]

# convert the sample sheet info into a pandas dataframe
sample_sheet = pd.read_csv("sample_sheet.csv")

# list of all of the samples from the sheet
samples = list(sample_sheet["sample"])

samples_and_ref = list(reference) + samples

# for nanopore data, aligner must be minimap2
sample_sheet["aligner"] = np.where(
    sample_sheet["platform"] == "nanopore", "minimap2", config["illumina"]["aligner"]
)

sample_sheet["trim"] = np.where(
    (
        pd.isnull(sample_sheet["read_1_adapter"].iloc[0])
        or pd.isnull(sample_sheet["read_2_adapter"].iloc[0])
    ),
    False,
    True,
)

chunks = list(range(1, config["freebayes_opts"]["chunks"] + 1))

if not "msa" in config["output"]:
    mask_repeats = config["mask_repeats"]
else:
    mask_repeats = False

# if not calling as groups, groups can be left blank in sample sheet
# if not, throw an error if any are left blank
try:
    group_set = set(list(sample_sheet["group"]))
    groups = {
        i: list(sample_sheet.loc[sample_sheet["group"] == i, "sample"])
        for i in list(group_set)
    }
except:
    if config["call_as_groups"]:
        raise ValueError(
            "call_as_groups option is enabled in config.yaml, but no groups are defined in the sample sheet"
        )


def is_paired(sample_name):
    sample = get_sample(sample_name)
    read_type = sample["read_type"].iloc[0].lower()
    if read_type == "paired":
        return True
    elif read_type == "single":
        return False
    else:
        raise ValueError("Read type must be either paired or single")


def get_read_source(sample):
    return sample["source"].iloc[0].lower()


def get_ref_basename(mask_repeats=mask_repeats):
    basename = "data/resources/"
    basename += config["ref_name"]
    if mask_repeats:
        basename += "_masked"
    return basename


# return reference file, name depends on whether repeats are being masked
def get_ref(w, fai=False):
    basename = get_ref_basename()
    basename += ".fa"
    if fai:
        basename += ".fai"
    return basename


def get_labels(reference=reference):
    with open(reference) as ref:
        pattern = re.compile(r"^>(\S+)")
        labels = []
        for line in ref:
            match = pattern.search(line)
            if match:
                labels.append(match.group(1))
    return labels


def get_fasta_to_mv(w):
    if w.sample_or_ref == config["ref_name"]:
        return f"data/resources/{config['ref_name']}.fa.split/{config['ref_name']}.part_{w.label}.fa"
    else:
        return f"data/consensus/{w.sample_or_ref}.fa.split/{w.sample_or_ref}.part_{w.label}.fasta"


def get_data_path(sample_name):
    sample = get_sample(sample_name)
    source = get_read_source(sample)
    location = sample["location"].iloc[0]
    if source == "novogene":
        return os.path.join(
            config["data_locations"][location],
            f"usftp21.novogene.com/01.RawData/",
            sample_name,
        )
    elif source == "nanopore":
        return os.path.join(
            config["data"][location], f"fastq_pass/{location}/", sample_name
        )
    else:
        raise ValueError(
            "One or more samples could not be found in the provided data location, but no valid alternate source was listed. The source in sample sheet must be either novogene or nanopore."
        )


def get_ids_for_sample(sample_name):
    path = get_data_path(sample_name)
    files = os.listdir(path)
    if is_paired(sample_name):
        regex = r"(?<=" + sample_name + r"_).*(?=_[12].fq.gz)"
    else:
        regex = r"(?<=" + sample_name + r"_).*(?=.fq.gz)"
    p = re.compile(regex)
    matches = [p.findall(i) for i in files]
    return [val for lst in matches for val in lst]


def get_alns_to_merge(w):
    return [
        os.path.join("data/alignments/", f"{w.sample}_{i}_sort.bam")
        for i in get_ids_for_sample(w.sample)
    ]


def get_reads(w):
    path = get_data_path(w.sample)
    if is_paired(w.sample):
        r1 = os.path.join(path, "{sample}_{iden}_1.fq.gz")
        r2 = os.path.join(path, "{sample}_{iden}_2.fq.gz")
        return [r1, r2]
    else:
        return os.path.join(path, "{sample}_{iden}.fq.gz")


def get_reads_to_trim(w):
    sample = get_sample(w.sample)
    if is_paired(w.sample):
        return [
            "data/reads/{sample}_{iden}_1.fq.gz",
            "data/reads/{sample}_{iden}_2.fq.gz",
        ]
    else:
        return "data/reads/{sample}_{iden}.fq.gz"


def get_sample(sample_name):
    return sample_sheet.loc[sample_sheet["sample"] == sample_name]


def get_adapter_seqs(w):
    sample = get_sample(w.sample)
    r1 = sample["read_1_adapter"]
    r2 = sample["read_2_adapter"]
    if r1.any():
        if r2.any():
            return f"-a {r1.iloc[0]} -A {r2.iloc[0]}"
        else:
            return f"-a {r1.iloc[0]}"
    elif r2.any():
        return f"-A {r2.iloc[0]}"


def get_fastqc_files(w):
    if is_paired(w.sample):
        return (
            expand(
                "data/qc/fastqc/{sample}_{read}_fastqc.zip", sample=samples, read=[1, 2]
            ),
        )
    else:
        return (expand("data/qc/fastqc/{sample}_fastqc.zip", sample=samples),)
    get_trim_files,


# get required names of aligned reads, determines which aligner is used
def get_aligned_reads(w):
    sample = get_sample(w.sample)
    if sample["aligner"].iloc[0].lower() == "bwa":
        return (f"data/alignments/{w.sample}_{w.iden}_bwa.bam",)
    elif sample["aligner"].iloc[0].lower() == "bowtie2":
        return (f"data/alignments/{w.sample}_{w.iden}_bt2.bam",)
    elif sample["aligner"].iloc[0].lower() == "minimap2":
        return (f"data/alignments/{w.sample}_{w.iden}_mm2.bam",)
    else:
        raise ValueError(
            "Invalid aligner in config.yaml: options are bwa, bowtie2 or minimap2"
        )


def get_alns_for_pileup(w, bai=False):
    if config["caller"] == "freebayes" or config["bcftools_opts"]["call_as_groups"]:
        return [f"data/alignments/{i}_dedup.bam" for i in groups[w.group]]
    else:
        return f"data/alignments/{w.sample_or_group}_dedup.bam"


# return either bams for all of the samples in the group or individual sample bams
def get_alns_in_group(w, groups=groups, bai=False):
    if config["bcftools_opts"]["call_as_groups"] or config["caller"] == "freebayes":
        return [f"data/alignments/{i}_dedup.bam" for i in groups[w.group]]
    else:
        return f"data/alignments/{w.sample}_dedup.bam"


# get name of variant caller from config file
def get_caller(w):
    return f"data/calls/{w.sample}_{config['caller']}_unprocessed.bcf"


def get_call_type(w):
    if config["caller"] == "freebayes" or config["bcftools_opts"]["call_as_groups"]:
        return f".tmp/group_call_{w.sample}_{config['caller']}.txt"
    else:
        return f".tmp/single_call_{w.sample}.txt"


# locate reads based on the location(s) listed in the sample sheet
def get_file_locations(w, read=None):
    sample = sample_sheet.loc[sample_sheet["sample"] == w.sample]
    base_dir = config["data_locations"][sample["location"].iloc[0]]
    if sample["read_type"].iloc[0] == "paired":
        return os.path.join(base_dir, w.sample + f"_{read}.fq.gz")
    elif sample["read_type"].iloc[0] == "single":
        return os.path.join(base_dir, w.sample + ".fq.gz")


# get the final bcf file, name depends on several options
def get_final_bcf(w, s=None, csi=False):
    if s == None:
        s = w.sample
    if not csi:
        return f"data/calls/{s}_norm_qflt.bcf"
    else:
        return f"data/calls/{s}_norm_qflt.bcf.csi"


# get final output for rule all
def get_final_output(w):
    out = []
    if "alignments" in config["output"]:
        [
            out.append(i)
            for i in expand("data/alignments/{sample}_dedup.bam.bai", sample=samples)
        ]
    if "calls" in config["output"]:
        [out.append(get_final_bcf(w, s=i, csi=True)) for i in samples]
    if "consensus" in config["output"]:
        [out.append(i) for i in expand("data/consensus/{sample}.fa", sample=samples)]
    if "msa" in config["output"]:
        [out.append(i) for i in get_msa_outputs()]
    if "tsvs" in config["output"]:
        out.append("data/tsvs/merged.tsv")
    out.append("data/qc/multiqc.html")
    return out


def get_fas_to_cat(w, ref=config["ref_name"], groups=groups):
    grp_and_ref = groups[f"{w.group}"] + [ref]
    return expand(
        "data/consensus/{sample_or_ref}_{{label}}_idfix.fa", sample_or_ref=grp_and_ref
    )


def get_msa_outputs(ref=reference, groups=groups, labels=get_labels()):
    # the reference should have the same names in the fasta as each of the consensus sequences
    outputs = []
    for i in groups:
        [outputs.append(f"data/msa/{j}_{i}.afa") for j in labels]
    return outputs


# return the group that a sample belongs to
def get_group_from_sample(w, groups=groups):
    for key, val in groups.items():
        if w.sample in val:
            return f"data/calls/{key}_{w.caller}_unprocessed.bcf"


# for ndj analysis, to merge progeny tsvs together
def get_progeny_tsvs(groups=groups):
    return ["data/tsvs/" + i + "_common.tsv" for i in groups["progeny"]]


# get quality cutoff for filtering
def get_qual_cutoff(groups=groups):
    return "-i'" + " & ".join(config["filtering"]) + "'"
    # caller = config["caller"]
    # cutoff_level = config["filtering"]["variant_quality_level"]
    # cutoff_val = config["filtering"]["variant_quality_cutoff_values"][caller][
    #     cutoff_level
    # ]
    # if caller == "freebayes":
    #     for key, val in groups.items():
    #         if w.sample in val:
    #             group = key
    #     n_samples = len(groups[group])
    #     cutoff_val = cutoff_val * n_samples
    # return cutoff_val


# get the format for transforming vcfs into tsvs
def get_query_format(w):
    cols = ["%SAMPLE", "%CHROM", "%POS", "%REF", "%ALT", "%QUAL", "%FILTER", "%GT"]
    if config["caller"] == "freebayes":
        add_cols = [
            "%DP",
            "%AF",
            "%AO",
            "%RO",
            "%TYPE",
        ]
    elif config["caller"] == "bcftools":
        add_cols = ["%DP", "%PL", "%DP4{0}", "%DP4{1}", "%DP4{2}", "%DP4{3}"]
    [cols.append(i) for i in add_cols]
    query = "\t".join(cols)
    return "[" + query + "]\n"


def get_reads_to_map(w, r=None):
    sample = sample_sheet.loc[sample_sheet["sample"] == w.sample]
    f = "data/reads/{sample}_{iden}"
    if sample["read_type"].iloc[0] == "paired":
        basename = [f"{f}_1", f"{f}_2"]
    else:
        basename = [f]
    if sample["trim"].any():
        final = [i + "_trimmed.fq.gz" for i in basename]
    else:
        final = [i + ".fq.gz" for i in basename]
    return final


# get ref index for bowtie2
def get_ref_bowtie2(w):
    return multiext(
        get_ref(w),
        ".1.bt2",
        ".2.bt2",
        ".3.bt2",
        ".4.bt2",
        ".rev.1.bt2",
        ".rev.2.bt2",
    )


# get ref index for bwa mem
def get_ref_bwa(w):
    return multiext(get_ref(w), ".amb", ".ann", ".bwt", ".pac", ".sa")


def get_ref_for_ref_parent(w, fai=False, mask_repeats=mask_repeats):
    if mask_repeats:
        base = f"data/resources/{config['ref_name']}_masked.fa"
    else:
        base = f"data/resources/{config['ref_name']}.fa"
    if fai:
        return base + ".fai"
    else:
        return base


# get ref index for minimap2
def get_ref_minimap2(w):
    return get_ref(w) + ".mmi"


# get ref file to generate freebayes regions
# no sample or group wildcard in this rule so relying on ref wildcard instead
def get_ref_to_generate_regions(w, fai=False, mask_repeats=mask_repeats):
    if mask_repeats:
        base = f"data/resources/{w.ref}_masked.fa"
    else:
        base = f"data/resources/{w.ref}.fa"
    if fai:
        base += ".fai"
    return base


# get a bed file for the reference being used to call variants for a particular sample
def get_region_from_sample(w, mask_repeats=mask_repeats):
    basename = get_ref(w)
    region_basename = re.sub(r"^(data/resources/)", r"\1regions/", basename)
    region_basename = re.sub(r"(_masked)?(\.fa)", r"", region_basename)
    return region_basename + f".{w.label}.region.{w.i}.bed"


# determines which chroms to use for variant calling
def get_regions_to_call(w):
    return ",".join(config["chroms"])


# get all of the tsv files to merge into one
def get_tsvs_to_merge(samples=samples):
    tsvs = [i for i in samples]
    return expand("data/tsvs/{sample}.tsv", sample=tsvs)
