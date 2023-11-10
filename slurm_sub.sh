#!/bin/bash

#SBATCH -p general
#SBATCH --mem=1G
#SBATCH -n 1
#SBATCH -t 5-
#SBATCH --mail-type=end
#SBATCH --mail-user=cannecar@ad.unc.edu

# shopt -s extglob
# module load anaconda
# conda create --name seq_pipeline python=3.11 snakemake
# module load bowtie2
# module load bwa
# module load samtools/1.18
# module load freebayes/1.1.0-54
# module load repeatmasker/4.1.5

# cd scripts || exit
# bash cat_lanes.sh
# cd .. || exit
# if [[ ! -f "input/reads/w1118_1.fq.gz" ]]; then
#   cd input || exit
#   prefetch SRR2044312
#   fasterq-dump SRR2044312
#   for i in SRR2044312*.f?(ast)q; do
#     bgzip $i
#   done
#   for i in [12]; do
#     mv SRR2044312_"$i".f?(ast)q.gz reads/w1118_"$i".fq.gz
#   done
#   cd .. || exit
# fi
# cp -n /proj/sekellab/projects/OregonR_Hawley_isogenized/oregonR.bam output/aligments/oregonr.bam
# cp -rn /proj/sekellab/projects/ndj_progeny/usftp21.novogene.com/01.RawData/ndj_+([0-9]) input/
module load repeatmasker
module load samtools
module load sambamba
snakemake --profile=slurm --use-conda
