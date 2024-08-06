#!/bin/bash

#SBATCH -p general
#SBATCH --mem=1G
#SBATCH -n 1
#SBATCH -t 5-
#SBATCH --mail-type=end
#SBATCH --mail-user=cannecar@ad.unc.edu

module purge
mamba activate ndj

export SNAKEMAKE_OUTPUT_CACHE=/proj/sekellab/.snakemake-cache/
snakemake --profile=slurm