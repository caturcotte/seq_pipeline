#!/bin/bash

#SBATCH -p general
#SBATCH --mem=1G
#SBATCH -n 1
#SBATCH -t 5-
#SBATCH --mail-type=end
#SBATCH --mail-user=cannecar@ad.unc.edu

module purge
module load seqkit
module load repeatmasker
module load samtools
module load sambamba
module load muscle
module load freebayes

snakemake --profile=slurm --use-conda
