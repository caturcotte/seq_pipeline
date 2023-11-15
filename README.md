# DNA sequencing pipeline

This pipeline is designed to take reads from Illumina or Nanopore data and generate any of the following:
- Alignments
- Variant calls
- Consensus sequences from variant calls
- TSVs to easily import variant calls into R for analysis

## Instructions
1. Log into longleaf and clone this repository to a location suitable for analysis (e.g. work space)
```
ssh onyen@longleaf.unc.edu
cd /work/users/o/n/onyen/
git clone https://github.com/caturcotte/seq_pipeline.git
```
### Editing the sample sheet and config files
2. In a separate terminal window (or exit your SSH session in the current terminal window by typing `exit`), cd to a directory you can easily find in your file browser. Download the sample sheet to this directory, open it in Excel, and edit it with your sample information.
- Adapters are optional. If left blank, nothing will be trimmed.
- `location` is just a nickname for a file path specified later in `config.yaml`.
```
cd ~/Downloads # or another directory that you can locate easily in Excel
sftp onyen@longleaf.unc.edu
cd /work/users/o/n/onyen/seq_pipeline/
get sample_sheet.csv
# Edit the sample sheet in Excel and save it, it should be in your Downloads or wherever you cd'd to initially
put sample_sheet.csv
```

3. Exit out of the SFTP connection and go back to your SSH connection to Longleaf. Enter the seq_pipeline directory created earlier. Edit the configuration file `config.yaml` based on the instructions in the file.
- For file locations: if the location of a sample called A in `sample_sheet.csv` is listed as `example`, and `example` is defined in `file_locations` in `config.yaml` as `/work/users/o/n/onyen/example/`, snakemake will look for the file `/work/users/o/n/onyen/example/A.fq.gz` (if single) or `/work/users/o/n/onyen/example/A_1.fq.gz` and `/work/users/o/n/onyen/example/A_2.fq.gz` (if paired) to start the pipeline.
```
cd seq_pipeline/
nano config.yaml # or vim config.yaml, or whatever command line text editor you like
```

4. Open `slurm/config.yaml` in a text editor and replace `YOUR EMAIL HERE` with your ONYEN email. Also enter your email into `--mail-user` in `slurm_sub.sh`.

### Installing mamba and activating the environment
   
6. Run `miniforge_installer.sh` - this will also create a conda environment to run snakemake called `seq_env`
```
cd seq_pipeline
bash miniforge_installer.sh
```

7. Add conda and mamba to your PATH.
```
source "~/conda/etc/profile.d/conda.sh"
source "~/conda/etc/profile.d/mamba.sh"
conda init
mamba init
source ~/.bashrc
```

8. Activate the `seq_env` environment created earlier.
```
mamba activate seq_env
```
### Running the pipeline

9. Perform a dry run of the pipeline with `snakemake -n`. It should spit out a list of jobs that it plans to run. If it fails, it's likely you did not input your samples correctly into the sample sheet.


10. If the dry run worked, submit the job to longleaf with `sbatch slurm_sub.sh`. You should be notified by email if the job completed or failed, but you can also check the status of your jobs on Longleaf with `squeue -u onyen`.
   
## Expected output
|**Output**|**File Location**|
|:--------:|:----------------|
|alignments|`data/alignments/{sample}.bam(.bai)`|
|calls|`data/calls/{sample}.bcf(.csi)`|
|consensus|`data/seqs/{sample}.fa`|
|tsvs|`data/tsvs/{sample}.tsv`, `merged.tsv` (all samples)|
