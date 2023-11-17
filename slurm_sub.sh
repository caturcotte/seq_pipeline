#!/bin/bash

#SBATCH -p general
#SBATCH --mem=1G
#SBATCH -n 1
#SBATCH -t 5-
#SBATCH --mail-type=end
#SBATCH --mail-user=cannecar@ad.unc.edu

module load repeatmasker
module load samtools
module load sambamba
snakemake --profile=slurm --use-conda
