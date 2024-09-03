# return sample row given sample name
def get_sample(sample_name, df):
    return df.loc[df["sample"] == sample_name]


# return split fields from sample name def get_sample_fields(self):
def get_fields(sample_name):
    fields = sample_name.split("-")
    return dict(zip(["condition", "sample_type", "sample_num"], fields))


# check if reads are illumina or nanopore
def is_short_read(sample_name, df):
    sample = get_sample(sample_name, df)
    platform = sample['platform'].iloc[0]
    if platform == "illumina":
        return True
    elif platform == "nanopore":
        return False
    else:
        raise ValueError(
            f"Invalid platform specified for sample {sample} in sample_sheet.csv. Platform must be either Illumina or Nanopore."
        )


def has_barcode(sample_name, df):
    sample = get_sample(sample_name, df)
    return sample["barcode"].notnull().all()


def get_barcode(sample_name, df):
    sample = get_sample(sample_name, df)
    return f"barcode{int(sample['barcode'].iloc[0]):02}"


# search each file in a directory for a regex match
def regex_over_dir(path, regex):
    p = re.compile(regex)
    matches = []
    for i in path.iterdir():
        match = p.search(str(i))
        if match is not None:
            matches.append(match)
    return matches


# return the path at which a sample's reads are located
def get_data_path(sample_name, df, config):
    sample = get_sample(sample_name, df)
    location = sample["location"].iloc[0]
    if location not in config["data_locations"]:
        raise ValueError(f"Sample location {location} not found in config file.")
    filepath = Path(config["data_locations"][location])
    if not filepath.exists():
        raise OSError(f"File directory {str(filepath)} not found.")
    # dirs = [f"usftp21.novogene.com/01.RawData/{sample_name}", sample_name]
    dirs=[sample_name]
    if has_barcode(sample_name, df):
        barcode = get_barcode(sample_name, df)
        # dirs.append(f"fastq_pass/{barcode:02}")
        date=pd.to_datetime(sample['date'].iloc[0]).strftime("%Y%m%d")
        d = f"{date}/"
    else:
        d = sample_name
    path = filepath.joinpath(d)
    if path.exists():
        return path
    else:
        raise OSError(f"Directories containing reads not found in {location}.")


# find all fastqs associated with a sample
def find_all_fastqs(sample_name, df, config):
    sample = get_sample(sample_name, df)
    path = get_data_path(sample_name, df, config)
    if has_barcode(sample_name, df): 
        regex = r"(\w*\.f(?:ast)?q(?:\.gz)?$)"
        # barcode = get_barcode(sample)
        # regex = rf"(FA[A-Z][0-9]{{5}}_pass_{barcode}_\w{{8}}_\w{{8}}_[0-9]*\.f(?:ast)?q(?:\.gz)?)"
    else:
        regex = rf"(.*{sample_name}_[a-zA-Z0-9_-]*(?:_[12])?\.f(?:ast)?q(?:\.gz)?)"
        # regex = rf"(.*{sample_name}_[a-zA-Z0-9_-]*(?:_[12])?\.f(?:ast)?q(?:\.gz)?)""
    matches = regex_over_dir(path, regex)
    if not matches:
        raise OSError(f"No reads found in {str(path)} for sample {sample_name}.")
    return [i[0] for i in matches]


def get_fastqs_with_gz(sample_name, df, config):
    fastqs = find_all_fastqs(sample_name, df, config)
    gz_matches = []
    for i in fastqs:
        if i.endswith("q"):
            gz = "".join([i, ".gz"])
        else:
            gz = i
        gz_matches.append(gz)
    return gz_matches


def get_fastq(sample_name, id, read, df, config):
    files = get_fastqs_with_gz(sample_name, df, config)
    path = get_data_path(sample_name, df, config)
    if has_barcode(sample_name, df):
        # if all(str(i).startswith("fastq") for i in files):
            # regex = rf"(\w*_{id}\.f(?:ast)?q\.gz)"
        # elif all(str(i).startswith("FA") for i in files):
        barcode = get_barcode(sample_name, df)
            # regex = rf"(FA[A-Z][0-9]{{5}}_pass_{barcode}_\w{{8}}_\w{{8}}_{id}\.f(?:ast)?q\.gz)"
        regex = rf"{barcode}.fq.gz"
    else:
        regex = rf"({sample_name}_){id}_[12]\.f(?:ast)?q\.gz)"
    matches = []
    p = re.compile(regex)
    for i in files:
        match = p.search(i)
        if match is not None:
            matches.append(match)
    if len(matches) > 1:
        raise OSError(f"More than one FASTQ found for sample {sample_name}, id {id}, read {read}")
    if matches:
        return [path.joinpath(i[0]) for i in matches]
    elif len(files) > 1:
        raise OSError(f"More than one FASTQ found for sample {sample_name}, id {id}, read {read}")
    return [path.joinpath(i) for i in files]


def get_ids_for_sample(sample_name, df, config):
    files = find_all_fastqs(sample_name, df, config)
    if has_barcode(sample_name, df):
        regex = r"(?:\w*_)([0-9]{1,5})\.f(?:ast)?q(?:\.gz)?"
        if all(str(i).startswith("fastq") for i in files):
            regex = r"(?:\w*_)([0-9]{1,5})\.f(?:ast)?q(?:\.gz)?"
        elif all(str(i).startswith("FA") for i in files):
            barcode = get_barcode(sample_name, df)
            regex = rf"(?<=FA[A-Z][0-9]{{5}}_pass_{barcode}_\w{{8}}_\w{{8}}_)([0-9]*)(?=(?:\.f(?:ast)?q(?:\.gz)?))"
    else:
        regex = rf"(?<={sample_name}_)([a-zA-Z0-9_-]*)(?=_[12]\.f(?:ast)?q(?:\.gz)?)"
    matches = []
    p = re.compile(regex)
    for i in files:
        match = p.search(i)
        if match is not None:
            matches.append(match)
    if matches:
        return [i[1] for i in matches]
    else:
        return [sample_name]


# check if a sample has lane IDs or not
def has_lane_id(sample_name, df, config):
    files = find_all_fastqs(sample_name, df, config)
    ids = get_ids_for_sample(sample_name, df, config)
    return bool(ids)


# get sample field information to create sample names
def concat_sample_names(df):
    df['sample_num'] = round(df['sample_num'])
    sample_fields = list(zip(
        df["condition"],
        df["sample_type"].fillna(''),
        df["sample_num"].fillna(''),
    ))
    s_fields = []
    for i in sample_fields:
        s_fields.append([j for j in i if j])
    for i in s_fields:
        if len(i) > 2:
            i[2] = str(int(i[2])).zfill(3)
    sample_names = ['-'.join(i) for i in s_fields]
    df['sample'] = sample_names
    return df


def get_sample_names(df):
    return list(df['sample'])


def get_parent_names(df, config):
    ref_parent_name = config["ref_parent"]
    alt_parent_name = config["alt_parent"]
    try:
        ref_parent = df.loc[df["sample"] == ref_parent_name]
        alt_parent = df.loc[df["sample"] == alt_parent_name]
    except:
        raise OSError("One or more parents were not found in the sample sheet.")
    return (ref_parent_name, alt_parent_name)


def get_progeny_names(df, config):
    return [i for i in get_sample_names(df) if i not in get_parent_names(df, config)]


def get_sample_type(sample_name, df, config):
    if sample_name in get_progeny_names(df):
        return 'progeny'
    if sample_name == config['ref_parent']:
        return 'ref_parent'
    if sample_name == config['alt_parent']:
        return 'alt_parent'
    raise ValueError(f'Sample {sample_name} not found in list of samples.')

# list all lane IDs for one sample
def list_all_idens(df, config):
    prefixes = []
    for sample_name in get_sample_names(df):
        [prefixes.append(i) for i in get_ids_for_sample(sample_name, df, config)]
    return prefixes


# get reference file name/base name, non input function, name depends on whether repeats are masked
def get_ref_no_input(config, base, fai):
    name = ["data/resources/", config["ref_name"]]
    if config['mask_repeats']:
        name.append("_masked")
    if base:
        return "".join(name)
    name.append(".fa")
    if not fai:
        return "".join(name)
    name.append(".fai")
    return "".join(name)
