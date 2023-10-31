rule bcftools_call_single:
    input:
        alns="mapped/{sample}_sort_dedup.bam",
        aln_idxs="mapped/{sample}_sort_dedup.bam.bai",
        ref=get_ref,
        ref_idx=lambda w: get_ref(w, fai=True),
    output:
        "called/{sample}_bcftools_unprocessed.bcf",
    params:
        regions=lambda w: get_regions_to_call(w),
    shell:
        "bcftools call -cv -f {input.ref} {input.alns} -r {params.regions} -o {output}"
