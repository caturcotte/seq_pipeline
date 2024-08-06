#!/bin/bash

set -euo pipefail

module load samtools
module load python/3.12.2
source .venv/bin/activate

basedir="/work/users/c/a/cannecar/parentals/seq_pipeline/data/calls/"
destdir="/work/users/c/a/cannecar/vcf2tiger/"
now="" &&
# echo "$now Submitting jobs to SLURM to convert VCFs to TSVs."
# for i in "${basedir}"*_bcftools_raw.bcf
# do
#   tempfile=${destdir}${i##*/}
#   destfile=${tempfile%_bcftools_raw.bcf}_snps.tsv
#   if [ -f "${destfile}" ]
#   then
#     echo "$now ${destfile##/} already exists. Skipping..."
#     continue
#   else
#     echo "$now Sumitting ""${i#"${basedir}"}"" for processing..."
#     sbatch -p general -n 1 -N 1 -t 3-0 --mem 50G --wait --wrap="./vcf_script.sh ${i} ${destfile}" > /dev/null &
#   fi
# done
#
# echo "$now Finished submitting jobs."
# echo "$now Waiting for jobs to finish..."
# wait &&
#
# echo "$now Finished converting VCFs to TSVs."
# echo "$now Moving files to correct destinations..."
# for i in {01..16} 
# do
#   mv "${destdir}/ndj_${i}_snps.tsv" "${destdir}/samples/WT-G0-0${i}.tsv" &
# done
#
# for i in oregonr_combined_snps.tsv w1118_combined_snps.tsv
# do
#   mv ${i} ${destdir}/references/${i%_combined_snps.tsv}_snps.tsv &
# done
# wait &&
#
# echo "$now Finished moving files."
echo "$now Converting TSVs to parquet..."

for file in samples/*.tsv references/*.tsv
do
    sbatch -p general --mem 24G -n 1 -t 5- --wait --wrap="python3 convert_tsv_to_parquet.py ${file}" > /dev/null &
done
wait &&

echo "$now Done!"
