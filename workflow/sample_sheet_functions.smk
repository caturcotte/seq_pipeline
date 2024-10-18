# return sample row given sample name
import sysconfig


def detect_platform():
    platform = sysconfig.get_platform()
    if platform.endswith("arm64"):
        architecture = "arm64"
    elif platform.endswith("universal2") or platform.endswith("x86_64"):
        architecture = "x64"
    if platform.startswith("linux"):
        os = "linux"
    elif platform.startswith("macos"):
        os = "macos"
    else:
        raise OSError("This pipeline can only run on Linux or MacOS systems.")
    return os, architecture


def get_pod5s(folders, config, df):
    base_dir = Path(config["base_dir"])
    pod5s = [str(i) for i in base_dir.glob("*/pod5*/*.pod5")]
    return pod5s


def get_sample(sample_name, df):
    return df.loc[df["sample"] == sample_name]


def concat_sample_names(df):
    df['sample_num_str'] = df['sample_num'].astype(str).str.zfill(3)
    df['sample'] = df[['genotype', 'sample_num_str']].agg('_'.join, axis=1)
    return df


# return split fields from sample name def get_sample_fields(self):
def get_fields(sample_name):
    fields = sample_name.split("-")
    return dict(zip(["condition", "sample_type", "sample_num"], fields))


# check if reads are illumina or nanopore
def is_short_read(sample, df):
    if sample in df['sample']:
        return False
    else:
        return True
    # sample = get_sample(sample_name, df)
    # platform = sample['platform'].iloc[0]
    # if platform == "illumina":
    #     return True
    # elif platform == "nanopore":
    #     return False
    # else:
    #     raise ValueError(
    #         f"Invalid platform specified for sample {sample} in sample_sheet.csv. Platform must be either Illumina or Nanopore."
    #     )


# def has_barcode(sample_name, df):
#     sample = get_sample(sample_name, df)
#     return sample["barcode"].notnull().all()


# def get_barcode(sample_name, df):
#     sample = get_sample(sample_name, df)
#     return f"barcode{int(sample['barcode'].iloc[0]):02}"


def get_parent_names(config):
    return {
        "reference": config["reference_genotype"]["name"],
        "variant": config["variant_genotype"]["name"],
        "paternal": config["paternal_genotype"]["name"],
    }


def get_all_sample_ids(df):
    return list(sample_sheet["sample_id"].unique())


def get_sample_id_from_folder(df, folder):
    df_flt = df.loc[df["ont_folder"] == folder]
    sample_id = list(df_flt["sample_id"].unique())
    if len(sample_id) > 1:
        samples_str = ", ".join(sample_id)
        raise OSError(f"{folder} found in multiple sample directories: {samples_str}")
    if sample_id is None:
        raise OSError(f"{folder} not found.")
    return sample_id[0]


def get_progeny_dict(df):
    progeny_dict = {}
    for i in list(df["ont_folder"].unique()):
        progeny_dict[i] = list(df.loc[df["ont_folder"] == i, "sample"])
    return progeny_dict


def get_progeny_names(df):
    return list(sample_sheet["sample"].unique())
