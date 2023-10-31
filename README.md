# DNA sequencing pipeline

This pipeline is designed to take reads from Illumina or Nanopore data and generate any of the following:
- Alignments
- Variant calls
- Consensus sequences from variant calls
- TSVs to easily import variant calls into R for analysis

## Instructions
1. Clone the repository to a suitable location for analysis. If working on longleaf:
```
cd /work/users/o/n/onyen/
git clone https://github.com/caturcotte/seq_pipeline.git
cd seq_pipeline/
```
2. Load conda and create the environment (you should only need to do this once):
```
module load anaconda
conda create -f workflow/envs/environment.yaml
```

3. Open `sample_sheet.xlsx` in Excel or a text editor and follow the instructions for inputting your sample information into the spreadsheet. Make sure you are entering data into the correct sheet (PAIRED_END for Illumina, SINGLE_END for Nanopore. Single end Illumina currently not supported.)
- If you want, you can copy the sample sheet to your host computer:
```
sftp onyen@longleaf.unc.edu
cd /work/users/o/n/onyen/seq_pipeline/
get sample_sheet.xlsx
# edit the file in Excel and save...
put sample_sheet.xlsx
```
- Or if you are on a Unix-based operating system, you can download the command line tool `sshfs` to temporarily mount the remote directory from longleaf to a location on your host computer:
```
mkdir seq_pipeline
sshfs onyen@longleaf.unc.edu:/work/users/o/n/onyen/seq_pipeline seq_pipeline/
```
- Then, locate `sample_sheet.xlsx` using your normal file browser, edit it in Excel and save.

Open `config.yaml` in a text editor and edit the options to suit your needs (instructions in file).

4. Open `slurm/config.yaml` in a text editor and replace `YOUR EMAIL HERE` with your ONYEN email. Also enter your email into `--mail-user` in `slurm_sub.sh`.

5. Perform a dry run of the pipeline with `snakemake -n`. If it fails, it's likely you did not input your samples correctly into the sample sheet.

6. If the dry run worked, submit the job to longleaf with `sbatch slurm_sub.sh`. You should be notified by email if the job completed or failed, but you can also check the status of your jobs on Longleaf with `squeue -u onyen`.
   
## Expected output
|**Output**|**File Location**|
|:--------:|:----------------|
|alignments|`mapped/{sample}.bam(.bai)`|
|calls|`called/{sample}.bcf(.csi)`|
|consensus|`seqs/{sample}.fa`|
|tsvs|`tsvs/{sample}.tsv`, `merged.tsv` (all samples)|
