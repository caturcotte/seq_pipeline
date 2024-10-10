#!/bin/bash

#SBATCH -p general
#SBATCH --mem=40G
#SBATCH -n 4
#SBATCH -t 5-
#SBATCH --mail-type=end
#SBATCH --mail-user=cannecar@ad.unc.edu

# module purge
# module load anaconda/2024.02
# conda activate ndj
# conda config --set channel_priority strict
# conda config --set solver libmamba

export SNAKEMAKE_OUTPUT_CACHE=/proj/sekellab/.snakemake-cache/
# snakemake --profile=slurm -p data/tiger/WT-F1-075/WT-F1-075.tiger_input.txt
 # snakemake -p "data/resources/chrom_lengths_tiger.txt" --profile=slurm
#  snakemake --profile=slurm -p data/calls/WT-F1-020_raw_clair3.bcf
snakemake --profile=slurm
 # snakemake --profile=slurm -p data/tiger/WT-F1-054/WT-F1-054_plot1.csv
