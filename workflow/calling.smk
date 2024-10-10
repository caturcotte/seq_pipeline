from pathlib import Path

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
        bcf="data/calls/{sample}_raw_bcftools.bcf",
    threads: 16
    envmodules:
        config['envmodules']['samtools']
    shell:
        "bcftools call -c -o {output.bcf} --threads {threads} {input.pileup}"


rule get_clair3_model:
    output:
        "workflow/scripts/rerio/clair3_models/r1041_e82_400bps_sup_v500_model"
    shell:
        "workflow/scripts/rerio/download_model.py rerio/clair3_models/r1041_e82_400bps_sup_v500_model"

# only call at candidate sites of parents to speed up
# min af of 0 makes it so that all calls including reference calls have information provided, because otherwise information won't be added when variant af is under threshold
rule clair3_call:
    input:
        alns="data/alignments/{sample}_markdup.bam",
        idxs="data/alignments/{sample}_markdup.bam.bai",
        # vcf="data/calls/parent_isec/unique_snps.vcf",
        # mdl="data/resources/r1041_e82_400bps_sup_v500/",
        ref=get_ref,
        ref_idx=lambda w: get_ref(w, fai=True)
    output:
        vcf="data/calls/{sample}/merge_output.vcf.gz",
    params:
        chromosomes=list(config['chromosomes'].keys()),
    conda:
        "envs/clair3.yaml"
    threads: 64
    shell:
        "run_clair3.sh -b {input.alns} -f {input.ref} -t {threads} -p ont -m /work/users/c/a/cannecar/co_analysis/seq_pipeline/data/resources/r1041_e82_400bps_sup_v500/ -o \"data/calls/{wildcards.sample}\" --ctg_name={params.chromosomes} --sample_name={wildcards.sample} --snp_min_af=0.0 --indel_min_af=0.0 --enable_phasing --use_whatshap_for_final_output_haplotagging --gvcf"

rule call_structural_variants:
    input:
        alns="data/alignments/{sample}_markdup.bam",
        idxs="data/alignments/{sample}_markdup.bam.bai",
        ref=get_ref,
        ref_idx=lambda w: get_ref(w, fai=True)
    output:
        vcf="data/calls/{sample}/{sample}_svs.vcf.gz",
    shell:
        "sniffles --input {input.alns} --reference {input.ref} -v {output}


rule switch_zygosity_based_on_sv:
    input:
        alns="data/alignments/{sample}_markdup.bam",
        idxs="data/alignments/{sample}_markdup.bam.bai",
        c3="data/calls/{sample}/merge_output.vcf.gz",
        sv="data/calls/{sample}/{sample}_svs.vcf.gz",
    output:
        "data/alignments/{sample}_clair3.vcf",
    threads: 8,
    shell:
        "pypy3 clair3.py SwitchZygosityBasedOnSVCalls --bam_fn {input.alns} --clair3_vcf_input {input.c3} --sv_vcf_input {input.sv} --vcf_output {output} --threads {threads}"
rule whatshap_phase:
    input:
        bam="data/alignments/{sample}_markdup.bam",
        idx="data/alignments/{sample}_markdup.bam.bai",
        vcf="data/calls/{sample}_raw_clair3.vcf",
        ref=get_ref,
        ref_idx=lambda w: get_ref(w, fai=True)
    output:
        "data/calls/{sample}_raw_phased.vcf"
    conda:
        "envs/whatshap.yaml"
    shell:
        "whatshap phase -o {output} --reference={input.ref} {input.vcf} {input.bam}"


rule unzip_and_mv_clair3:
    input:
        "data/calls/{sample}/merge_output.vcf.gz",
    output:
        "data/calls/{sample}_raw_clair3.vcf",
    envmodules:
        config['envmodules']['samtools']
    shell:
        "bgzip -d {input} -o {output}"
        

rule convert_phased_vcf:
    input:
        "data/calls/{sample}_raw_phased.vcf"
    output:
        "data/calls/{sample}_raw_phased.bcf"
    envmodules:
        config['envmodules']['samtools']
    shell:
        "bcftools view -Ob -o {output} {input}"