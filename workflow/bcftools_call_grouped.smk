rule bcftools_call_grouped:
    input:
        alns=get_alns_for_calling,
        aln_idxs=lambda w: get_alns_for_calling(w, bai=True),
        ref=get_ref,
        ref_idx=lambda w: get_ref(w, fai=True),
    output:
        "called/{group}_bcftools_unprocessed.bcf",
    shell:
        "bcftools call -cv -f {input.ref} {input.alns} -o {output}"

