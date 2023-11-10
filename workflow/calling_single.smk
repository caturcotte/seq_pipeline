rule bcftools_mpileup_single:
    input:
        alns="mapped/{sample}_sort_dedup.bam",
        ref=get_ref,
        ref_idx=lambda w: get_ref(w, fai=True),
    output:
        "called/{sample}_bcftools_pileup.vcf.gz",
    params:
        regions=lambda w: get_regions_to_call(w),
        extra="--max-depth 200 --min-BQ 15"
    threads: 16
    conda:
        "envs/bcftools.yaml"
    shell:
        "bcftools mpileup -f {input.ref} -r {params.regions} {params.extra} -o {output} --threads {threads} {input.alns}"

rule bcftools_call_single:
    input:
        "called/{sample}_bcftools_pileup.vcf.gz",
    output:
        bcf="called/{sample}_bcftools_unprocessed.bcf",
        call_type=temp('.tmp/bcftools_single.txt')
    params:
        regions=lambda w: get_regions_to_call(w),
    threads: 4
    conda:
        "envs/bcftools.yaml"
    shell:
        "bcftools call -cv -r {params.regions} -o {output.bcf} --threads {threads} {input} && touch {output.call_type}"

