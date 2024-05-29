import pandas as pd
import numpy as np
import re
from pathlib import Path


reference = config["reference"]

# convert the sample sheet info into a pandas dataframe
sample_sheet = pd.read_csv("sample_sheet.csv")

# list of all of the samples from the sheet
samples = list(sample_sheet["sample"])

samples_and_ref = list(reference) + samples

# for nanopore data, aligner must be minimap2
sample_sheet["aligner"] = np.where(
    sample_sheet["platform"] == "nanopore", "minimap2", config["alignment"]["aligner"]
)

sample_sheet["trim"] = np.where(
    (
        pd.isnull(sample_sheet["read_1_adapter"].iloc[0])
        or pd.isnull(sample_sheet["read_2_adapter"].iloc[0])
    ),
    False,
    True,
)

chunks = list(range(1, config["calling"]["freebayes_opts"]["chunks"] + 1))

if config["mask_repeats"] and "msa" in config["output"]:
    raise ValueError(
        "Masking repeats is not compatible with MSA output. Please remove MSA from required output or disable mask_repeats in config.yaml."
    )

# if not calling as groups, groups can be left blank in sample sheet
# if not, throw an error if any are left blank
try:
    group_set = set(list(sample_sheet["group"]))
    groups = {
        i: list(sample_sheet.loc[sample_sheet["group"] == i, "sample"])
        for i in list(group_set)
    }
except:
    # if config["call_as_groups"]:
    raise ValueError(
        " in the sample sheet."
    )


all_suffix_options = ["biallelic_", "snps_", "homozygous_", "flt_"]
suffix_constraints = []
for i in range(1, len(all_suffix_options)+1):
    suffix_constraints_i = combinations(all_suffix_options, i)
    suffix_constraints.append(suffix_constraints_i)
suffix_constraints = list(chain.from_iterable(suffix_constraints))
suffix_constraints = ["".join(list(i)) for i in suffix_constraints]
suffix_constraints = [i.rstrip("_") for i in list(suffix_constraints)]


def is_paired(sample_name):
    sample = get_sample(sample_name)
    read_type = sample["read_type"].iloc[0].lower()
    if read_type == "paired":
        return True
    elif read_type == "single":
        return False
    else:
        raise ValueError(
            f"Incorrect read type specified for sample {sample} in sample_sheet.csv. Read type must be either paired or single."
        )


def get_read_source(sample):
    return sample["source"].iloc[0].lower()


# return reference file, name depends on whether repeats are being masked
def get_ref(w=None, base=False, fai=False):
    name = ["data/resources/", config["ref_name"]]
    if config["mask_repeats"]:
        name.append("_masked")
    if base:
        return "".join(name)
    else:
        name.append(".fa")
        if fai:
            name.append(".fai")
        return "".join(name)


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
        return (
            "data/consensus/{sample_or_ref}.fa.split/{sample_or_ref}.part_{label}.fasta"
        )


def get_data_path(sample_name):
    sample = get_sample(sample_name)
    source = get_read_source(sample)
    location = sample["location"].iloc[0]
    filepath = Path(config["data_locations"][location])
    if source == "novogene":
        return filepath.joinpath(
            f"usftp21.novogene.com/01.RawData/",
            sample_name,
        )
    elif source == "nanopore":
        return filepath.joinpath(f"fastq_pass/{location}/", sample_name)
    else:
        raise ValueError(
            "One or more samples could not be found in the provided data location, but no valid alternate source was listed. The source in sample sheet must be either novogene, nanopore or other."
        )


def get_ids_for_sample(sample_name):
    path = get_data_path(sample_name)
    if not check_for_lane_id(sample_name):
        return ["nolane"]
    else:
        if is_paired(sample_name):
            regex = r"(?<=" + sample_name + r"_).*(?=_[12].fq.gz)"
        else:
            regex = r"(?<=" + sample_name + r"_).*(?=.fq.gz)"
        p = re.compile(regex)
        matches = [p.findall(str(i)) for i in path.iterdir()]
        return [val for lst in matches for val in lst]


def get_alns_to_merge(w):
    return [
        str(Path("data/alignments/").joinpath(f"{w.sample}_{i}_sort.bam"))
        for i in get_ids_for_sample(w.sample)
    ]


def check_for_lane_id(sample_name):
    sample = get_sample(sample_name)
    if sample["lane_ids"].iloc[0]:
        return True
    else:
        return False


def get_reads(w):
    path = get_data_path(w.sample)
    sample = get_sample(w.sample)
    if check_for_lane_id(w.sample):
        prefix = "{sample}_{iden}"
    else:
        prefix = "{sample}"
    if is_paired(w.sample):
        r1 = path.joinpath(f"{prefix}_1.fq.gz")
        r2 = path.joinpath(f"{prefix}_2.fq.gz")
        return [str(i) for i in [r1, r2]]
    else:
        return str(path.joinpath(f"{prefix}.fq.gz"))


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
            "Invalid aligner in config.yaml: options are bwa, bowtie2 or minimap2."
        )


def get_alns_for_pileup(w, bai=False):
      return [f"data/alignments/{i}_markdup.bam" for i in groups[w.group]]


# return either bams for all of the samples in the group or individual sample bams
def get_alns_in_group(w, groups=groups, bai=False):
    # if config["bcftools_opts"]["call_as_groups"] or config["caller"] == "freebayes":
    return [f"data/alignments/{i}_markdup.bam" for i in groups[w.group]]
    # else:
    #     return f"data/alignments/{w.sample}_dedup.bam"


def get_regions():
    if (
        "chromosomes" in config["calling"].keys()
        or "regions" in config["calling"].keys()
    ):
        return "-R data/resources/call_regions.bed"
    else:
        return ""


# get name of variant caller from config file
def get_caller(w):
    return f"data/calls/{w.sample}_{config['caller']}_unprocessed.bcf"


def group_or_single(w):
    # if config["calling"]["bcftools_opts"]["call_as_groups"]:
    return "data/calls/{sample_or_group}_bcftools_group_pileup.bcf"
    # else:
    #     return "data/calls/{sample_or_group}_bcftools_single_pileup.bcf"


# locate reads based on the location(s) listed in the sample sheet
def get_file_locations(w, read=None):
    sample = sample_sheet.loc[sample_sheet["sample"] == w.sample]
    base_dir = Path(config["data_locations"][sample["location"].iloc[0]])
    if sample["read_type"].iloc[0] == "paired":
        return base_dir.join(base_dir, w.sample + f"_{read}.fq.gz")
    elif sample["read_type"].iloc[0] == "single":
        return base_dir.join(base_dir, w.sample + ".fq.gz")


def get_filtering_criteria(w):
    all_options = {
        "min_depth": "FORMAT/DP < ",
        "min_gq": "FORMAT/GQ < ",
        "min_qd": "INFO/QD < ",
    }
    flt_strs = []
    for key, val in config["filtering"]["all"].items():
        if key in all_options:
            flt_str = "".join([all_options[key], str(val)])
            flt_strs.append(flt_str)
    if w.group in config["filtering"]:
        for key, val in config["filtering"][w.group].items():
            if key in all_options:
                flt_str = "".join([all_options[key], str(val)])
                flt_strs.append(flt_str)
    return " & ".join(flt_strs)


def get_bcf_suffix(
    sample_or_group,
    caller=config['calling']['caller'],
    suffix_list=all_suffix_options,
):
    all_options = {
        "biallelic_only": "biallelic",
        "snps_only": "snps",
        "homozygous_only": "homozygous",
        "min_depth": "flt",
        "min_gq": "flt",
        "min_qd": "flt",
        "raw_options": "flt",
    }
    suffixes = []
    suffix_list = [i.rstrip("_") for i in suffix_list]
    print(suffix_list)
    for opt in all_options:
        if opt in config["filtering"]["all"]:
            if config["filtering"]["all"][opt]:
                suffixes.append(all_options[opt])
        if sample_or_group in config["filtering"]:
            if opt in config["filtering"][sample_or_group]:
                if config["filtering"][sample_or_group][opt]:
                    suffixes.append(all_options[opt])

    # remove duplicate values but keep ordered
    suffixes = sorted(list(set(suffixes)), key=lambda x: suffix_list.index(x))
    suffix_str = "_".join(suffixes)
    return suffix_str


# get the final bcf file, name depends on several options
def get_final_bcf(sample_or_group, csi=False, caller=config['calling']['caller']):
    suffix_str = get_bcf_suffix(sample_or_group)
    if not csi:
        return f"data/calls/{sample_or_group}_{caller}_{suffix_str}.bcf"
    else:
        return f"data/calls/{sample_or_group}_{caller}_{suffix_str}.bcf.csi"

def get_call_file(w, addon=None, csi=False):
    calls = "".join(["data/calls/{group}_", f"{config['calling']['caller']}_", "{suffix}.bcf"])
    if csi:
        calls = "".join([calls, ".csi"])
    return calls

def get_bcf_wildcard(w, csi=False):
    return get_final_bcf(w.group, csi)

# get final output for rule all
def get_final_output(w, samples=samples, groups=groups):
    out = []
    samples_or_groups = []
    # for i in config["analysis"]["comparisons"]:
    #     for j in config["analysis"]["comparisons"][i]["datasets"]:
    #         samples_or_groups.append(j)
    for i in config["filtering"]:
        if i != "all":
           samples_or_groups.append(i) 
    samples_and_groups = list(set(list(groups.keys()) + samples_or_groups))
    if "alignments" in config["output"]:
        [
            out.append(i)
            for i in expand("data/alignments/{sample}_markdup.bam.bai", sample=samples)
        ]
    if "calls" in config["output"]:
        [out.append(get_final_bcf(i, csi=True)) for i in samples_and_groups]
    if "consensus" in config["output"]:
        [out.append(i) for i in expand("data/consensus/{sample}.fa", sample=samples)]
    if "msa" in config["output"]:
        [out.append(i) for i in get_msa_outputs()]
    if "tsvs" in config["output"]:
        [out.append(i) for i in expand("data/tsvs/{group}.tsv", group=groups.keys())]
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


def get_group_name(w, groups=groups):
    for key, val in groups.items():
        if w.sample in val:
            return key


# return the group that a sample belongs to
def get_group_from_sample(w, groups=groups):
    return get_final_bcf(get_group_name(w))
    # if config["filtering"]["snps_only"]:
    #     return f"data/calls/{get_group_name(w)}_{w.caller}_snps.bcf"
    # else:
    #     return f"data/calls/{get_group_name(w)}_{w.caller}.bcf"


# for ndj analysis, to merge progeny tsvs together
def get_progeny_tsvs(groups=groups):
    return ["data/tsvs/" + i + "_common.tsv" for i in groups["progeny"]]


# get quality cutoff for filtering
def get_qual_cutoff(groups=groups):
    filters = " & ".join(config["filtering"])
    return "".join(["-i'", filters, "'"])
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


def get_comparisons(w, basename=False, csi=False):
    comp = config["analysis"]["comparisons"][{w.comparison}]
    if basename:
        return comp
    comp_files = []
    for i in comp:
        suffix = get_suffix_str(i)
        comp_file = "".join(["data/calls/", i, suffix, ".bcf"])
        comp_files.append(comp_file)
    if not csi:
        return comp_files
    return ["".join([i, ".csi"]) for i in comp_files]

def get_suffix_for_comparison(w):
    return get_suffix_str(w.f1)

def get_final_bams(w, bai=False):
    sample = get_sample(w.sample)
    file = ["data/alignments/{sample}"]
    if sample["PCR"]:
        file.append("_markdup")
    file.append(".bam")
    if bai:
        file.append(".bai")
    return "".join(file)


# get the format for transforming vcfs into tsvs
def get_query_format(w, caller=config['calling']['caller']):
    cols = ["%SAMPLE", "%CHROM", "%POS", "%REF", "%ALT", "%QUAL", "%FILTER", "%GT"]
    if caller in ["freebayes", "both"]:
        add_cols = [
            "%DP",
            "%AF",
            "%AO",
            "%RO",
            "%TYPE",
        ]
    elif caller in ["bcftools", "both"]:
        add_cols = ["%DP", "%PL", "%DP4{0}", "%DP4{1}", "%DP4{2}", "%DP4{3}"]
    [cols.append(i) for i in add_cols]
    query = "\t".join(cols)
    return "[" + query + "]\n"


def get_reads_to_map(w, r=None):
    sample = get_sample(w.sample)
    if check_for_lane_id(w.sample):
        f = "data/reads/{sample}_{iden}"
    else:
        f = "data/reads/{sample}_nolane"
    if sample["read_type"].iloc[0] == "paired":
        basename = [f"{f}_1", f"{f}_2"]
    else:
        basename = [f]
    if sample["trim"].any():
        final = ["".join([i, "_trimmed.fq.gz"]) for i in basename]
    else:
        final = ["".join([i, ".fq.gz"]) for i in basename]
    return final


# get ref index for bowtie2
def get_ref_bowtie2(w):
    return multiext(
        get_ref(w), ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"
    )


# get ref index for bwa mem
def get_ref_bwa(w):
    return multiext(get_ref(w), ".amb", ".ann", ".bwt", ".pac", ".sa")


def get_ref_for_ref_parent(w, fai=False):
    if config["mask_repeats"]:
        base = f"data/resources/{config['ref_name']}_masked.fa"
    else:
        base = f"data/resources/{config['ref_name']}.fa"
    if fai:
        return "".join([base, ".fai"])
    else:
        return base


# get ref index for minimap2
def get_ref_minimap2(w):
    return "".join([get_ref(w), ".mmi"])


# get ref file to generate freebayes regions
# no sample or group wildcard in this rule so relying on ref wildcard instead
def get_ref_to_generate_regions(w, fai=False):
    if config["mask_repeats"]:
        base = "data/resources/{ref}_masked.fa"
    else:
        base = "data/resources/{ref}.fa"
    if fai:
        base = "".join([base, ".fai"])
    return base


# get a bed file for the reference being used to call variants for a particular sample
def get_region_from_sample(w):
    basename = get_ref(w)
    region_basename = re.sub(r"^(data/resources/)", r"\1regions/", basename)
    if config["mask_repeats"]:
        region_basename = re.sub(r"(_masked)?(\.fa)", r"", region_basename)
    return "".join([region_basename, ".{label}.region.{i}.bed"])


# def get_ref_freebayes(fai=False):
#     if config["calling"]["chromosomes"] or config["calling"]["regions"]:
#         if not fai:
#             return "data/resources/call_regions.fa"
#         else:
#             return "data/resources/call_regions.fa.fai"


# get all of the tsv files to merge into one
# def get_tsvs_to_merge(samples=samples):
#     tsvs = [i for i in samples]
#     return expand("data/tsvs/{sample}.tsv", sample=tsvs)


def get_tsvs_to_merge(w):
    comparisons = {}
    for key, val in config["comparisons"]:
        if val["comparison"] == "unique":
            comparisons[key] == f"{val['datasets'][0]}_not_in_{val['datasets'][1]}"
        elif val["comparison"] == "common":
            comparisons[key] == f"{val['datasets'][0]}_in_{val['datasets'][1]}"
    return expand(
        "data/tsvs/{comparison}_{condition}.tsv",
        comparison=comparisons.keys(),
        condition=comparisons.values(),
    )
