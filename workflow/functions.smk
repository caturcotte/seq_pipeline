import pandas as pd
import re


reference = config["reference"]

# chunks to use for parallelizing freebayes
chunks = list(range(1, config["freebayes_opts"]["nchunks"] + 1))

# chroms to call variants for
chroms = config["chroms"]

# if doing ndj analysis, the ref for snp calling is one of the parentals
# the snps of that ref are determined by calling variants against dm6
# then a consensus sequence is built from that
if config["ndj_analysis"]:
    ref_names = [
        config["ref_name"],
        config["ndj_analysis_opts"]["parents"]["ref_parent"],
    ]
else:
    ref_names = [config["ref_name"]]

# paired end or single end is determined by seq type, not including
# Illumina single end as an option since we don't use it
if config["sequencer"] == "nanopore":
    seq_type = "SINGLE_END"
else:
    seq_type = "PAIRED_END"

# convert the sample sheet info into a pandas dataframe
sample_sheet = pd.read_excel("sample_sheet.xlsx", sheet_name=seq_type)

# list of all of the samples from the sheet
samples = list(sample_sheet["sample"])

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


# get required names of aligned reads, determines which aligner is used
def get_aligned_reads(w):
    if config["aligner"].lower() == "bwa":
        return (f"mapped/{w.sample}_bwa.sam",)
    elif config["aligner"].lower() == "bowtie2":
        return (f"mapped/{w.sample}_bt2.sam",)
    elif config["aligner"].lower() == "minimap2":
        return (f"mapped/{w.sample}_mm2.sam",)
    else:
        raise ValueError(
            "Invalid aligner in config.yaml: options are bwa, bowtie2 or minimap2"
        )


# return either bams for all of the samples in the group or individual sample bams
def get_alns_for_calling(w, groups=groups, bai=False):
    if config["bcftools_opts"]["call_as_groups"]:
        return [f"mapped/{i}_sort_dedup.bam" for i in groups[w.group]]
    else:
        return f"mapped/{w.sample}_sort_dedup.bam"


# get name of variant caller from config file
def get_caller(w):
    return f"called/{w.sample}_{config['caller']}_unprocessed.bcf"


# locate reads based on the location(s) listed in the sample sheet
def get_file_locations(w, end=None):
    sample = sample_sheet.loc[sample_sheet["sample"] == w.sample]
    if config["sequencer"] != "nanopore":
        if end == "r1":
            return sample["path_to_reads_1"]
        elif end == "r2":
            return sample["path_to_reads_2"]
    else:
        return sample["path_to_reads"]


# get the final bcf file, name depends on several options
def get_final_bcf(w, s=None, csi=False):
    if s == None:
        s = w.sample
    prefix = [f"called/{s}_norm"]
    if config["filtering"]["variant_quality"] != "off":
        prefix.append("qflt")
    if config["ndj_analysis"]:
        if s in config["ndj_analysis_opts"]["parents"].values():
            prefix.append("het")
    if not csi:
        return "_".join(prefix) + ".bcf"
    else:
        return "_".join(prefix) + ".bcf" + ".csi"


# get final output for rule all
def get_final_output(w):
    out = []
    if config["output"] == "tsvs":
        out.append("tsvs/merged.tsv")
        if config["ndj_analysis"]:
            out.append("tsvs/progeny_common.tsv")
        return out
    elif config["output"] == "consensus":
        return expand("seqs/{sample}.fa", sample=samples)
    elif config["output"] == "calls":
        vcfs = []
        for i in samples:
            vcfs.append(get_final_bcf(w, s=i, csi=True))
        return vcfs
    elif config["output"] == "alignments":
        return expand("mapped/{sample}_sort_dedup.bam.bai", sample=samples)


# return the group that a sample belongs to
def get_group_from_sample(w, groups=groups):
    for key, val in groups.items():
        if w.sample in val:
            return f"called/{key}/{w.caller}_unprocessed.bcf"


# optical duplicate distance for removing optical duplicates from alignments
# values are from samtools markdup docs
def get_opt_dup_distance(wildcards, config=config):
    if config["sequencer"].lower() == "hiseq":
        return "1000"
    elif config["sequencer"].lower() == "novaseq":
        return "2500"
    else:
        raise ValueError(
            "Invalid sequencer in config.yaml, options are HiSeq, NovaSeq or Nanopore"
        )


# for ndj analysis, to merge progeny tsvs together
def get_progeny_tsvs(w, groups=groups):
    return ["tsvs/" + i + "_common.tsv" for i in groups["progeny"]]


# get quality cutoff for filtering
def get_qual_cutoff(w):
    cutoff_level = config["filtering"]["variant_quality_level"]
    caller = config["caller"]
    return config["filtering"]["variant_quality_cutoff_values"][caller][cutoff_level]


# get the format for transforming vcfs into tsvs
def get_query_format(w):
    cols = ["%SAMPLE", "%CHROM", "%POS", "%REF", "%ALT", "%QUAL", "%FILTER", "%GT"]
    if config["caller"] == "freebayes":
        add_cols = [
            "%FORMAT/DP",
            "%AF_TOTAL",
            "%AO_TOTAL" "%DP_TOTAL",
            "%RO_TOTAL",
            "%TYPE",
            "%DP",
            "%RO",
            "%AO",
        ]
    elif config["caller"] == "bcftools":
        add_cols = ["%DP", "%PL", "%DP4{0}", "%DP4{1}", "%DP4{2}", "%DP4{3}"]
    [cols.append(i) for i in add_cols]
    return r"\t".join(cols) + r"\n"


# return reference file, name depends on whether repeats are being masked
def get_ref(w, fai=False):
    basename = "resources/"
    if not config["ndj_analysis"]:
        basename += config["ref_name"]
    else:
        try:
            s = w.group
        except:
            s = w.sample
        if s == config["ndj_analysis_opts"]["parents"]["ref_parent"]:
            basename += "dm6"
        else:
            basename += config["ndj_analysis_opts"]["parents"]["ref_parent"]
    final = [basename, "fa"]
    # final.append("fa")
    if config["mask_repeats"]:
        final.append("masked")
    if fai:
        final.append("fai")
    return ".".join(final)


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


# get ref index for minimap2
def get_ref_minimap2(w):
    return get_ref(w) + ".mmi"


# get ref file to generate freebayes regions
# no sample or group wildcard in this rule so relying on ref wildcard instead
def get_ref_to_generate_regions(w, fai=False):
    if config["mask_repeats"]:
        base = f"resources/{w.ref}.fa.masked"
    else:
        base = f"resources/{w.ref}.fa"
    if fai:
        base += ".fai"
    return base


# get a bed file for the reference being used to call variants for a particular sample
def get_region_from_sample(w):
    basename = get_ref(w)
    region_basename = re.sub(r"(^resources/)", r"\1regions/", basename)
    if config["mask_repeats"]:
        region_basename = re.sub(r"(.fa)(.masked)*", r"", region_basename)
    return region_basename + f".{w.chrom}.region.{w.i}.bed"


# determines which chroms to use for variant calling
def get_regions_to_call(w):
    return ",".join(config["chroms"])


# get all of the tsv files to merge into one
def get_tsvs_to_merge(w, samples=samples):
    if config["ndj_analysis"]:
        tsvs = [
            i
            for i in samples
            if i != config["ndj_analysis_opts"]["parents"]["ref_parent"]
        ]
    else:
        tsvs = [i for i in samples]
    return expand("tsvs/{sample}.tsv", sample=tsvs)

