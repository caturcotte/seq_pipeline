# DNA sequencing pipeline

This pipeline is designed to take reads from Illumina or Nanopore data and generate any of the following:
- Alignments
- Variant calls
- Consensus sequences from variant calls
- TSVs to easily import variant calls into R for analysis

## Instructions
1. Log into longleaf and download and run the conda installer miniforge:
 https://github.com/conda-forge/miniforge/releases

 ```
 ssh onyen@longleaf.unc.edu
 wget  https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge-pypy3-Linux-x86_64.sh
 sh Miniforge-pypy3-Linux-x86_64.sh
 ```

2. Create and activate the mamba environment needed to run snakemake.
 ```
 mamba create -n seq_pipeline python=3.11 snakemake=7.32.4
 mamba activate seq_pipeline
 ```
3. Clone the repository to a suitable location for analysis. If working on longleaf:
```
cd /work/users/o/n/onyen/
git clone https://github.com/caturcotte/seq_pipeline.git
cd seq_pipeline/
```

4. Open `sample_sheet.csv` in Excel or a text editor and put your sample information into the spreadsheet. Adapters are optional, if left blank nothing will be trimmed.

- If you want, you can copy the sample sheet to your host computer:
```
sftp onyen@longleaf.unc.edu
cd /work/users/o/n/onyen/seq_pipeline/
get sample_sheet.csv
# edit the file in Excel and save...
put sample_sheet.csv
```
- Or if you are on a Unix-based operating system, you can download the command line tool `sshfs` to temporarily mount the remote directory from longleaf to a location on your host computer:
```
mkdir seq_pipeline
sshfs onyen@longleaf.unc.edu:/work/users/o/n/onyen/seq_pipeline seq_pipeline/
```
- Then, locate `sample_sheet.csv` using your normal file browser, edit it in Excel and save.

- For file locations, use a shorthand. You can then define this shorthand in config.yaml. For instance, if the location of a sample called A in sample_sheet.csv is listed as proj, and proj is defined in config.yaml as /proj/sekellab/projects/example_project/, snakemake will look for the file /proj/sekellab/projects/example_project/A.fq.gz (if single) or /proj/sekellab/example_project/A_1.fq.gz and A_2.fq.gz (if paired) to start the pipeline.

5. Open `config.yaml` in a text editor and edit the options to suit your needs (instructions in file).

6. Open `slurm/config.yaml` in a text editor and replace `YOUR EMAIL HERE` with your ONYEN email. Also enter your email into `--mail-user` in `slurm_sub.sh`.

7. Perform a dry run of the pipeline with `snakemake -n`. If it fails, it's likely you did not input your samples correctly into the sample sheet.

8. If the dry run worked, submit the job to longleaf with `sbatch slurm_sub.sh`. You should be notified by email if the job completed or failed, but you can also check the status of your jobs on Longleaf with `squeue -u onyen`.
   
## Expected output
|**Output**|**File Location**|
|:--------:|:----------------|
|alignments|`data/alignments/{sample}.bam(.bai)`|
|calls|`data/calls/{sample}.bcf(.csi)`|
|consensus|`data/seqs/{sample}.fa`|
|tsvs|`data/tsvs/{sample}.tsv`, `merged.tsv` (all samples)|
