class Sample:
    def __init__(self, name, sheet):
        self.name = name
        self.sheet = sheet


    # return sample row given sample name
    def get_sample(self):
        return self.sheet.df.loc[self.sheet.df["sample"] == self.name]


    def populate_attrs_from_sheet(self):
        sample = self.get_sample()
        self.condition = sample['condition']
        self.platform = sample['platform']
        self.location = sample['location']
        if self.name in self.sheet.get_progeny_names():
            self.generation = sample['generation']
            self.sample_num = sample['sample_num'] 
        if self.platform == 'nanopore':
            self.flow_cell = sample['flow_cell']
            self.date = sample['date']
            self.fly = sample['fly']


    # return split fields from sample name
    def get_sample_fields(self):
        fields = self.split("-")
        return dict(zip(["condition", "sample_type", "sample_num"], fields))


    # check if reads are illumina or nanopore
    def is_short_read(self):
        if platform == "illumina":
            return True
        elif platform == "nanopore":
            return False
        else:
            raise ValueError(
                f"Invalid platform specified for sample {sample} in sample_sheet.csv. Platform must be either Illumina or Nanopore."
            )


    # check if a sample has lane IDs or not
    def has_lane_id(self):
        sample = self.get_sample()
        files = self.find_all_fastqs()
        ids = self.get_ids_for_sample()
        return bool(ids)


    def has_barcode(self):
        sample = self.get_sample()
        return sample["barcode"].notnull().all()


    def get_barcode(sample_name):
        sample = get_sample(sample_name)
        return f"barcode{int(sample['barcode'].iloc[0]):02}"


    # return the path at which a sample's reads are located
    def get_data_path(sample_name):
        sample = get_sample(sample_name)
        location = sample["location"].iloc[0]
        if location not in config["data_locations"]:
            raise ValueError(f"Sample location {location} not found in config file.")
        filepath = Path(config["data_locations"][location])
        if not filepath.exists():
            raise OSError(f"File directory {str(filepath)} not found.")
        # dirs = [f"usftp21.novogene.com/01.RawData/{sample_name}", sample_name]
        dirs=[sample_name]
        if has_barcode(sample_name):
            barcode = get_barcode(sample_name)
            # dirs.append(f"fastq_pass/{barcode:02}")
            date=pd.to_datetime(sample['date']).strftime("%Y%m%d")
            d = f"{date}/{barcode}"
        else:
            d = sample_name
        path = filepath.joinpath(d)
        if path.exists():
            return path
        else:
            raise OSError(f"Directories containing reads not found in {location}.")


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
        if has_barcode(sample_name): 
            barcode = get_barcode(sample_name)
            regex = rf"(FA[A-Z][0-9]{{5}}_pass_{barcode}_\w{{8}}_\w{{8}}_[0-9]*\.f(?:ast)?q(?:\.gz)?)"
        else:
            regex = rf"(.*{sample_name}_[a-zA-Z0-9_-]*(?:_[12])?\.f(?:ast)?q(?:\.gz)?)"
        matches = regex_over_dir(path, regex)
        if not matches:
            raise OSError(f"No reads found in {str(path)} for sample {sample_name}.")
        return [i[0] for i in matches]


    def get_ids_for_sample(sample_name):
        files = find_all_fastqs(sample_name)
        if has_barcode(sample_name):
            barcode = get_barcode(sample_name)
            regex = rf"(?<=FA[A-Z][0-9]{{5}}_pass_{barcode}_\w{{8}}_\w{{8}}_))([0-9]*)(?=(?:\.f(?:ast)?q(?:\.gz)?))"
        else:
            regex = rf"(?<={sample_name}_))([a-zA-Z0-9_-]*)(?=_[12]\.f(?:ast)?q(?:\.gz)?)"
        matches = regex_over_dir(path, regex)
        if matches:
            return [i[0] for i in matches]
        else:
            return [sample_name]


class SampleSheet:
    def __init__(self, df, config):
        self.df = df
        self.config = config


    def add_sample_names(self):
        # get sample field information to create sample names
        self.df['sample_num'] = round(self.df['sample_num'])
        sample_fields = list(zip(
            self.df["condition"],
            self.df["sample_type"].fillna(''),
            self.df["sample_num"].fillna(''),
        ))
        s_fields = []
        for i in sample_fields:
            s_fields.append([j for j in i if j])
        for i in s_fields:
            if len(i) > 2:
                i[2] = str(int(i[2])).zfill(3)
        sample_names = ['-'.join(i) for i in s_fields]
        self.df['sample'] = sample_names


    def get_samples(self):
        samples = []
        for i in self.sample_sheet['sample']
        s = Sample(i, self)
        s.populate_attrs_from_sheet()
        samples.append(s)
        self.samples = samples


    def get_sample_names(self):
        return list(sample_sheet['sample'])


    def get_parent_names(self):
        self.ref_parent = self.config["ref_parent"]
        self.alt_parent = self.config["alt_parent"]
        try:
            ref_parent = self.sample_sheet.loc[sample_sheet["sample"] == self.ref_parent]
            alt_parent = self.sample_sheet.loc[sample_sheet["sample"] == self.alt_parent]
        except:
            raise OSError("One or more parents were not found in the sample sheet.")


    def get_progeny_names(self):
        self.progeny = [i for i in self.get_sample_names() if i not in self.get_parent_names()]

    def get_sample(self, sample_name):
        return *[s for s in self.samples if s.name == sample_name]
    
    def check_property(self, sample_name, func):
        sample = self.get_sample(sample_name)
        return sample.func()
        

    # list all lane IDs for one sample
    def list_all_idens(sample_sheet):
        prefixes = []
        for sample_name in list(sample_sheet["sample"]):
            [prefixes.append(i) for i in get_ids_for_sample(sample_name)]
        return prefixes



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


    def has_barcode(sample_name):
        sample = self.get_sample(sample_name)
        return sample["barcode"].notnull().all()


    def get_barcode(self, sample_name):
        sample = self.get_sample(sample_name)
        return f"barcode{int(sample['barcode'].iloc[0]):02}"


    # return the path at which a sample's reads are located
    def get_data_path(self, sample_name):
        sample = self.get_sample(sample_name)
        location = sample.location
        if location not in config["data_locations"]:
            raise ValueError(f"Sample location {location} not found in config file.")
        filepath = Path(config["data_locations"][location])
        if not filepath.exists():
            raise OSError(f"File directory {str(filepath)} not found.")
        # dirs = [f"usftp21.novogene.com/01.RawData/{sample_name}", sample_name]
        dirs=[sample_name]
        if has_barcode(sample_name):
            barcode = get_barcode(sample_name)
            # dirs.append(f"fastq_pass/{barcode:02}")
            date=pd.to_datetime(sample['date']).strftime("%Y%m%d")
            d = f"{date}/{barcode}"
        else:
            d = sample_name
        path = filepath.joinpath(d)
        if path.exists():
            return path
        else:
            raise OSError(f"Directories containing reads not found in {location}.")


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
        if has_barcode(sample_name): 
            barcode = get_barcode(sample_name)
            regex = rf"(FA[A-Z][0-9]{{5}}_pass_{barcode}_\w{{8}}_\w{{8}}_[0-9]*\.f(?:ast)?q(?:\.gz)?)"
        else:
            regex = rf"(.*{sample_name}_[a-zA-Z0-9_-]*(?:_[12])?\.f(?:ast)?q(?:\.gz)?)"
        matches = regex_over_dir(path, regex)
        if not matches:
            raise OSError(f"No reads found in {str(path)} for sample {sample_name}.")
        return [i[0] for i in matches]


    def get_ids_for_sample(sample_name):
        files = find_all_fastqs(sample_name)
        if has_barcode(sample_name):
            barcode = get_barcode(sample_name)
            regex = rf"(?<=FA[A-Z][0-9]{{5}}_pass_{barcode}_\w{{8}}_\w{{8}}_))([0-9]*)(?=(?:\.f(?:ast)?q(?:\.gz)?))"
        else:
            regex = rf"(?<={sample_name}_))([a-zA-Z0-9_-]*)(?=_[12]\.f(?:ast)?q(?:\.gz)?)"
        matches = regex_over_dir(path, regex)
        if matches:
            return [i[0] for i in matches]
        else:
            return [sample_name]

def make_sample_class(sheet_file, config):
    sample_sheet = pd.read_csv(sheet_file)
    validate(sample_sheet, "schemas/sample.schema.yaml")
    all_samples = Samples(sample_sheet, config)
    all_samples.add_sample_names()
    all_samples.get_samples()
    all_samples.get_parents()
    all_samples.get_progeny()
    return all_samples