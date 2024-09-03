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
snakemake --version

export SNAKEMAKE_OUTPUT_CACHE=/proj/sekellab/.snakemake-cache/
snakemake --profile=slurm --notemp
