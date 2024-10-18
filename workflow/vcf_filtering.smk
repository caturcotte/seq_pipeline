# for some reason you can't just select for snps and the reference so instead you just have to exclude everything else
# this weirdly doesn't get rid of all indels so I remove them when rearranging the parquets
rule select_snps:
    input:
        "data/calls/{sample}_raw.bcf"
        # bcf="data/calls/{sample}_raw_{caller}.bcf",
        # csi="data/calls/{sample}_raw_{caller}.bcf.csi",
    output:
        "data/calls/{sample}.bcf",
    threads: 4
    envmodules:
        config['envmodules']['samtools']
    shell:
        "bcftools view -V mnps,bnd,indels,other -Ob -o {output} {input[0]}"
# rule get_unique_snps:
#     input:
#         ref=f"data/calls/{config['ref_parent']}.bcf",
#         ref_idx=f"data/calls/{config['ref_parent']}.bcf.csi",
#         alt=f"data/calls/{config['alt_parent']}.bcf",
#         alt_idx=f"data/calls/{config['alt_parent']}.bcf.csi"
#     output:
#         expand("data/calls/parent_isec/000{i}.vcf", i=[0, 1, 2, 3])
#     shell:
#         "bcftools isec -p data/calls/parent_isec/ {input.ref} {input.alt}"


# rule compress_unique_snps:
#     input:
#         "data/calls/parent_isec/000{i}.vcf",
#     output:
#         "data/calls/parent_isec/000{i}.vcf.gz"
#     envmodules:
#         config['envmodules']['samtools']
#     shell:
#         "bgzip {input}"

# rule index_unique_snps:
#     input:
#         "data/calls/parent_isec/000{i}.vcf.gz"
#     output:
#         "data/calls/parent_isec/000{i}.vcf.gz.csi"
#     envmodules:
#         config['envmodules']['samtools']
#     shell:
#         "bcftools index {input}"

# rule merge_unique_snps:
#     input:
#         vcf1="data/calls/parent_isec/0000.vcf.gz",
#         vcf2="data/calls/parent_isec/0001.vcf.gz",
#         csi1="data/calls/parent_isec/0000.vcf.gz.csi",
#         csi2="data/calls/parent_isec/0001.vcf.gz.csi"
#     output:
#         "data/calls/parent_isec/unique_snps.vcf"
#     envmodules:
#         config['envmodules']['samtools']
#     shell:
#         "bcftools merge {input.vcf1} {input.vcf2} -o {output}"

