#!/usr/bin/env bash

module load python/3.12.2
pip3 install pyarrow

basedir="data/calls/"
destdir="data/tsvs/"
# mkdir -p data/tsvs/references
# mkdir -p data/tsvs/references
# echo "Submitting jobs to SLURM to convert VCFs to TSVs."
# for i in "${basedir}"*_bcftools_raw.bcf
# do
#   tempfile=${destdir}${i##*/}
#   destfile=${tempfile%_bcftools_raw.bcf}_snps.tsv
#   if [ -f "${destfile}" ]
#   then
#     echo "${destfile##/} already exists. Skipping..."
#     continue
#   else
#     echo "Sumitting ""${i#"${basedir}"}"" for processing..."
#     sbatch -p general -n 1 -N 1 -t 3-0 --mem 50G --wait --wrap="./vcf_script.sh ${i} ${destfile}" > /dev/null &
#   fi
# done

# echo "Finished submitting jobs."
# echo "Waiting for jobs to finish..."
# wait &&

# echo "Finished converting VCFs to TSVs."
# echo "Moving files to correct destinations..."
# for i in {01..16} 
# do
#   # mv "${destdir}ndj_${i}_snps.tsv" "${destdir}WT-G0-0${i}.tsv"
#   sbatch -p general -n 1 -N 1 -t 3-0 --mem 20G --wait --wrap="python3 convert_tsv_to_parquet.py ${destdir}WT-G0-0${i}.tsv" &
# done
for i in oregonr_combined w1118_combined 
do
  # mv "${destdir}ndj_${i}_snps.tsv" "${destdir}WT-G0-0${i}.tsv"
  sbatch -p general -n 1 -N 1 -t 3-0 --mem 20G --wait --wrap="python3 convert_tsv_to_parquet.py ${destdir}${i}_snps.tsv" &
done
