
> [!important] 
> When a line in a code block starts with a `$`, it means that that line should be entered into the terminal prompt. **Do not** include the `$` when typing/copying the command into the prompt.
> 
> I will also embed anything that should be filled in in `<these brackets>`. Don't include the brackets when you fill in your own information.
> 
> **I HIGHLY recommend learning some [basic terminal skills](https://www.freecodecamp.org/news/linux-command-line-tutorial/) before attempting this tutorial.** I will show step by step instructions for how to run the pipeline, but will not include a tutorial for basic terminal navigation skills. This is just one guide, but there are many guides out there you can use as a reference.

>[!note]
>This guide is designed for use with Unix systems (Mac, Linux). Windows users may need to use [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/install) for some commands to work.
## Downloading data

Before starting with the pipeline, check out the guide to [[Uploading sequencing data to Longleaf|uploading your data to Longleaf]].
## Cloning the pipeline

If you haven't already, go to the terminal and connect to Longleaf via SSH:

```
$ ssh <onyen>@longleaf.unc.edu
```

Go to your work directory on Longleaf at `/work/users/o/n/onyen/` (where 'o' and 'n' are the first and second letters of your ONYEN, respectively).

```
$ cd /work/users/<o>/<n>/<onyen>/
```

Clone the pipeline from [here](https://github.com/caturcotte/seq_pipeline/tree/main) (click the green button that says Code, then copy the link in the drop down), then go to the directory of the repository.

```
$ git clone https://github.com/caturcotte/seq_pipeline.git
$ cd seq_pipeline/
```

The directory tree of the pipeline should look like this:

```
$ tree
.
├── config.yaml
├── miniforge_installer.sh
├── README.md
├── sample_sheet.csv
├── slurm
│   └── config.yaml
├── slurm_sub.sh
├── Snakefile
└── workflow
    ├── alignment.smk
    ├── calling.smk
    ├── envs
    │   ├── bcftools.yaml
    │   ├── bowtie2.yaml
    │   ├── bwa.yaml
    │   ├── concat_vcfs.yaml
    │   ├── environment.yaml
    │   ├── freebayes.yaml
    │   ├── minimap2.yaml
    │   ├── repeatmasker.yaml
    │   └── samtools.yaml
    ├── functions.smk
    ├── misc.smk
    ├── qc.smk
    └── vcf_filtering.smk
```

Now that the pipeline has been copied to Longleaf, we can configure it.

## Configuring the pipeline

### Editing `config.yaml`

>[!note]
There are two files called config.yaml - one in seq_pipeline/ and one in seq_pipeline/slurm/. Here I am referring to the one in seq_pipeline/.

Open config.yaml in a text editor:
```
$ nano config.yaml
```

You will see that there are numerous options for configuring the pipeline for your needs. There are comments in the file explaining the purpose of all of the options. Some of the more notable ones:
- **`ref_name` and `ref_file`:** The name and location of the reference genome [[FASTA]], respectively. The default is the dm6 release from UCSC, which is already stored on Longleaf. The pipeline will [symlink](https://en.wikipedia.org/wiki/Symbolic_link) (make a shortcut) to your ref_file at `data/resources/ref_name.fa`. Note that the original files are not altered in the pipeline, so there isn't a need to make a copy.
- **`data_locations`:** These are shorthands for specific locations where your data is located. See the [[Running the DNA sequencing analysis pipeline#Pointing the pipeline to your read files|Pointing the pipeline to your read files]] section for more detail.

### Adding your samples into sample_sheet.csv

In a new terminal, move to a **directory you can easily access via your file browser**. In this directory, open up an [[SFTP]] connection to [[Longleaf]] and copy `sample_sheet.csv` to your local computer.
```
$ cd ~/Downloads
$ sftp <onyen>@longleaf.unc.edu
$ cd /work/users/<o>/<n>/<onyen>/seq_pipeline/
$ get sample_sheet.csv
```

In this case, the file was downloaded to the Downloads folder.

>[!tip]
>If you open your SFTP connection in the wrong location, you can use `lls` (local list directory) and `lcd` (local change directory) to navigate to a better location without exiting the connection.

Edit sample_sheet.csv in Excel. Each of your samples will have its own row in the sheet. 

>[!note]
The adapter fields (Novogene and Nanopore) and barcode field (for Novogene only) can be left blank. If left blank, adapters/barcodes will not be trimmed for that sample. Often reads come demultiplexed, meaning reads containing adapter/barcodes have already been trimmed because the barcode sequences were used to sort the files into the folders for their respective samples.

**For Novogene samples** the sample names must match the original names you gave the samples when you sent them to Novogene. 

>[!tip]
>
You can check your sample names with `ls /path/to/reads/usftp21.novogene.com/01.RawData/` - all of the subdirectories in this directory will correspond to the names you gave your samples originally. 

**For Nanopore samples** the sample name can be whatever you want so long as it does not include spaces, but the barcode column should contain the number of the barcode used for that sample.

Save `sample_sheet.csv` as a csv file under the same name and directory as before (rewrite the original file). Then re-upload the file to Longleaf:

```
$ put sample_sheet.csv
```

You can now exit the SFTP connection (type `exit`) and go back to your original SSH connection (if you saved the other terminal window, otherwise connect to Longleaf and change directory to the `seq_pipeline` folder again).

### Pointing the pipeline to your read files

>[!note]
Multiple locations can be specified in data_locations. Each location should be on a new line in the file.

The way to point the pipeline to your files differs slightly between Novogene and Nanopore data.

Here I show examples for both. In both cases we will assume the `data_locations` entry in `config.yaml` looks like this, with an entry called `example_location`:

```
data_locations:
	example_location: /absolute/path/to/files/
```
#### For Novogene

If, in `sample_sheet.csv`, the location of sample A is set to `example_location`:

| sample | group  | platform | source   | read_type | location |
| ------ | ------ | -------- | -------- | --------- | -------- |
| A      | group1 | illumina | novogene | paired    | example_location  |

The pipeline will look for these files for A:
`/absolute/path/to/files/usftp21.novogene.com/01.RawData/A/A*_1.fq.gz` (read 1) and `/absolute/path/to/files/usftp21.novogene.com/01.RawData/A/A*_2.fq.gz` (read 2).

>[!note]
`*` is a [wildcard](https://en.wikipedia.org/wiki/Wildcard_character) meaning "one or more occurances of any character." For example, `A*_1.fq.gz` will match all read 1 files for sample A.

#### For Nanopore
If, in `sample_sheet.csv`, the location of sample A is set to `example_location`:

| sample | group  | platform | source   | read_type | location | barcode |
| ------ | ------ | -------- | -------- | --------- | -------- |
| A      | group1 | nanopore | nanopore | single    | example_location  | 1 |

The pipeline will look for the files for A in `/absolute/path/to/files/fastq_pass/barcode01/*.fq.gz`.

### Modifying slurm submission scripts with your email

Open slurm_sub.sh:
```
$ nano slurm_sub.sh
```

Edit the mail-user line, replacing `YOUR-EMAIL-HERE` with your UNC email.
Save and exit the file with `Ctrl+O` and `Ctrl+X`.

Do the same for `slurm/config.yaml`.

## Running the pipeline

We're ready to run the pipeline! Make sure you're in `seq_pipeline`, then use the following command to run the pipeline:

```
$ sbatch slurm_sub.sh
```

The output should say `Submitted batch job JOB-NUMBER.`

The pipeline is now running. While it's running, you can [[Monitor the status of your jobs|monitor the status of your jobs]] and read up on [[How the pipeline works|how the pipeline works]].

## Viewing your output files

All of the output files should be in the `data/` directory:

| output type | location         | recommended tool to view |
| ----------- | ---------------- | ------------------------ |
| alignments  | data/alignments/ | IGV                      |
| calls       | data/calls/      | text editor/less/bcftools              |
| consensus   | data/consensus/  | Snapgene alignment/IGV                         |
| tsvs            |  data/tsvs                | R                         |

You can use IGV on Longleaf directly so long as you use the `-X` flag when establishing your SSH connection:

```
ssh -X <onyen>@longleaf.unc.edu
```

>[!note]
>To use the `-X` flag on Mac you'll need to download [XQuartz](https://www.xquartz.org/).


## Checking the quality of your samples

The pipeline will automatically run several QC metrics for your samples and compile them into one report. To view the report, it's easiest to download the html content first to your home computer.

Exit your SSH connection and connect to Longleaf via SFTP. Then, navigate to the `data/qc` folder within `seq_pipeline` and download the multiqc files:
```
$ sftp <onyen>@longleaf.unc.edu
$ cd /work/users/<o>/<n>/<onyen>/seq_pipeline/data/qc/
$ get -r multiqc*

# files will download...

$ exit
```

Now, open `multiqc.html` in Firefox:
```
$ firefox multiqc.html
```

There are various metrics in this report that you can look at. A good one to start with is the [[FASTQC]] report, which tells you the quality of your reads. You can also see the average read depth for each of your samples, and [[Transition to Transversion Ratio]] to determine the quality of the calls. For *Drosophila* WGS, the Ts/Tv ratio should be around 2 [1]. Too low indicates a high false positive rate of calling, and too high indicates a biased callset (see [this post](**https://bioinformatics.stackexchange.com/questions/4362/why-ti-tv-ratio**))

See [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035531572-Evaluating-the-quality-of-a-germline-short-variant-callset) for more information on variant calling metrics.

[1] Keightley PD, Trivedi U, Thomson M, Oliver F, Kumar S, Blaxter ML. Analysis of the genome sequences of three Drosophila melanogaster spontaneous mutation accumulation lines. Genome Res. 2009 Jul;19(7):1195-201. doi: 10.1101/gr.091231.109. Epub 2009 May 13. PMID: 19439516; PMCID: PMC2704435.