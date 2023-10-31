# DNA sequencing pipeline

This pipeline is designed to take reads from Illumina or Nanopore data and generate any of the following:
- Alignments
- Variant calls
- Consensus sequences from variant calls
- TSVs to easily import variant calls into R for analysis

## Instructions
Clone the repository to a suitable location for analysis. If working on longleaf:
```
cd /work/users/o/n/onyen/
git clone https://github.com/caturcotte/seq_pipeline.git
cd seq_pipeline/
```
Open `sample_sheet.xlsx` in Excel or a text editor and follow the instructions for inputting your sample information into the spreadsheet. Make sure you are entering data into the correct sheet (PAIRED_END for Illumina, SINGLE_END for Nanopore. Single end Illumina currently not supported.)

Open `config.yaml` in a text editor and edit the options to suit your needs (instructions in file).

Open `slurm/config.yaml` in a text editor and replace `YOUR EMAIL HERE` with your ONYEN email. Also enter your email into `--mail-user` in `slurm_sub.sh`.

Perform a dry run of the pipeline with `snakemake -n`. If it fails, it's likely you did not input your samples correctly into the sample sheet.

If the dry run worked, submit the job to longleaf with `sbatch slurm_sub.sh`. You should be notified by email if the job completed or failed, but you can also check the status of your jobs on Longleaf with `squeue -u onyen`.
