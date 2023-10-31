rule bcftools_call_single:
    input:
        alns="mapped/{sample}_sort_dedup.bam",
        aln_idxs="mapped/{sample}_sort_dedup.bam.bai",
        ref=get_ref,
        ref_idx=lambda w: get_ref(w, fai=True),
    output:
        "called/{sample}_bcftools_unprocessed.bcf",
    shell:
        "bcftools call -cv -f {input.ref} {input.alns} -o {output}"
