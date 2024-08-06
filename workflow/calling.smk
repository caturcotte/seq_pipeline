rule bcftools_mpileup:
    input:
        alns="data/alignments/{sample}_markdup.bam",
        idxs="data/alignments/{sample}_markdup.bam.bai",
        ref=get_ref,
        ref_idx=lambda w: get_ref(w, fai=True),
    output:
        bcf=temp("data/calls/{sample}_pileup.bcf"),
    threads: 32
    envmodules:
        config['envmodules']['samtools']
    shell:
        "bcftools mpileup -f {input.ref} -a AD,DP -Ou -o {output.bcf} --threads {threads} {input.alns}"


rule bcftools_call:
    input:
        pileup="data/calls/{sample}_pileup.bcf",
    output:
        bcf="data/calls/{sample}_raw.bcf",
    threads: 16
    envmodules:
        config['envmodules']['samtools']
    shell:
        "bcftools call -c -o {output.bcf} --threads {threads} {input.pileup}"
