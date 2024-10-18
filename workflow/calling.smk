from pathlib import Path


rule bcftools_mpileup:
    input:
        aln=multiext("data/alignments/{parent}", ".bam", ".bam.bai"),
        ref=multiext(
            "data/resources/genome",
            ".fa",
            ".fa.fai",
        ),
    output:
        bcf=temp("data/calls/{parent}_pileup.bcf"),
    threads: 32
    envmodules:
        config["envmodules"]["samtools"],
    shell:
        "bcftools mpileup -f {input.ref[0]} -a AD,DP -Ou -o {output.bcf} --threads {threads} {input.aln[0]}"


rule bcftools_call:
    input:
        pileup="data/calls/{parent}_pileup.bcf",
    output:
        bcf="data/calls/{parent}_raw.bcf",
    threads: 16
    envmodules:
        config["envmodules"]["samtools"],
    shell:
        "bcftools call -c -o {output.bcf} --threads {threads} {input.pileup}"


rule clair3_call:
    input:
        aln=multiext("data/alignments/{progeny}", ".bam", ".bam.bai"),
        ref=multiext(
            "data/resources/genome",
            ".fa",
            ".fa.fai",
        ),
        mdl="data/resources/rerio/clair3_models/r1041_e82_400bps_sup_v500/pileup.index",
    output:
        "data/calls/{progeny}/merge_output.vcf.gz",
    params:
        chromosomes=list(config["chromosomes"].keys()),
        out_dir=lambda w, output: str(Path(output[0]).parents[0])
    conda:
        "envs/clair3.yaml"
    threads: 64
    shell:
        'run_clair3.sh -b {input.aln[0]} -f {input.ref[0]} -t {threads} -p ont -m /work/users/c/a/cannecar/co_analysis/seq_pipeline/data/resources/r1041_e82_400bps_sup_v500/ -o {params.out_dir} --ctg_name={params.chromosomes} --sample_name={wildcards.progeny} --snp_min_af=0.0 --indel_min_af=0.0 --enable_phasing --gvcf'


# rule call_structural_variants:
#     input:
#         aln=multiext("data/alignments/{progeny}", ".bam", ".bam.bai"),
#         ref=multiext(
#             "data/resources/genome",
#             ".fa",
#             ".fa.fai",
#         ),
#     output:
#         vcf="data/calls/{progeny}_svs.vcf.gz",
#     shell:
#         "sniffles --input {input.aln[0]} --reference {input.ref[0]} -v {output}"


# rule switch_zygosity_based_on_sv:
#     input:
#         aln=multiext("data/alignments/{progeny}", ".bam", ".bam.bai"),
#         c3="data/calls/{progeny}/merge_output.vcf.gz",
#         sv="data/calls/{progeny}_svs.vcf.gz",
#     output:
#         "data/calls/{progeny}_raw.vcf",
#     threads: 8
#     shell:
#         "pypy3 clair3.py SwitchZygosityBasedOnSVCalls --bam_fn {input.aln[0]} --clair3_vcf_input {input.c3} --sv_vcf_input {input.sv} --vcf_output {output} --threads {threads}"


# rule whatshap_phase:
#     input:
#         bam="data/alignments/{sample}_markdup.bam",
#         idx="data/alignments/{sample}_markdup.bam.bai",
#         vcf="data/calls/{sample}_raw_clair3.vcf",
#         ref=get_ref,
#         ref_idx=lambda w: get_ref(w, fai=True)
#     output:
#         "data/calls/{sample}_raw_phased.vcf"
#     conda:
#         "envs/whatshap.yaml"
#     shell:
#         "whatshap phase -o {output} --reference={input.ref} {input.vcf} {input.bam}"


rule unzip_and_mv_clair3:
    input:
        "data/calls/{sample}/merge_output.vcf.gz",
    output:
        "data/calls/{sample}_raw.bcf",
    envmodules:
        config['envmodules']['samtools']
    shell:
        "bgzip -d {input} -c | bcftools view -Ob -o {output}"
# rule convert_phased_vcf:
#     input:
#         "data/calls/{sample}_raw_phased.vcf"
#     output:
#         "data/calls/{sample}_raw_phased.bcf"
#     envmodules:
#         config['envmodules']['samtools']
#     shell:
#         "bcftools view -Ob -o {output} {input}"
