dir="data/calls/${3}"
bcftools isec -p ${dir} ${1} ${2}
bcftools view "${dir}/0000.vcf" -Ob -o "${dir}/${4}_not_in_${5}.bcf"
bcftools view "${dir}/0001.vcf" -Ob -o "${dir}/${5}_not_in_${4}.bcf"
bcftools view "${dir}/0002.vcf" -Ob -o "${dir}/${4}_in_${5}.bcf"
bcftools view "${dir}/0003.vcf" -Ob -o "${dir}/${5}_in_${4}.bcf"
