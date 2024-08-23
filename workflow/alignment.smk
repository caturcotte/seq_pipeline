# ruleorder: merge_bams > mv_nolane_bams


rule align_bowtie2:
    input:
        reads=get_reads_to_map,
        ref=get_ref_bowtie2,
    output:
        temp("data/alignments/{sample}_{iden}_bt2.sam"),
    params:
        ref_basename=lambda w: get_ref(w, base=True),
    threads: 32
    conda:
        "envs/bowtie2.yaml"
    shell:
        "bowtie2 -x {params.ref_basename} -1 {input.reads[0]} -2 {input.reads[1]} -p {threads} --rg-id {wildcards.iden} --rg 'SM:{wildcards.sample}' > {output}"


rule align_minimap2:
    input:
        reads=get_reads_to_map,
        ref=get_ref_minimap2,
    output:
        temp("data/alignments/{sample}_{iden}_mm2.sam"),
    threads: 32
    envmodules:
        config['envmodules']['minimap2']
    shell:
        "minimap2 -ax map-ont {input.ref} -R '@RG\\tID:{wildcards.iden}\\tSM:{wildcards.sample}' -t {threads} -o {output} {input.reads}"


rule sam2bam:
    input:
        "data/alignments/{sample}_{iden}_{aligner}.sam",
    output:
        temp("data/alignments/{sample}_{iden}_{aligner}.bam"),
    envmodules:
        config['envmodules']['samtools']
    shell:
        "samtools view -1 -o {output} {input}"


rule fix_mate_pairs:
    input:
        get_aligned_reads,
    output:
        temp("data/alignments/{sample}_{iden}.bam"),
    envmodules:
        config['envmodules']['samtools']
    shell:
        "samtools fixmate -m -O bam,level=1 {input} {output}"


rule sort_bams:
    input:
        "data/alignments/{sample}_{iden}.bam",
    output:
        temp("data/alignments/{sample}_{iden}_sort.bam"),
    threads: 8
    envmodules:
        config['envmodules']['samtools']
    shell:
        "samtools sort {input} -l 1 -o {output} --threads {threads}"


# rule mv_nolane_bams:
#     input:
#         "data/alignments/{sample}_{sample}_sort.bam",
#     output:
#         temp("data/alignments/{sample}.bam"),
#     shell:
#         "mv {input} {output}"


rule merge_bams:
    input:
        get_alns_to_merge,
    output:
        temp("data/alignments/{sample}.bam"),
    threads: 4
    envmodules:
        config['envmodules']['sambamba']
    script:
        "scripts/mv_or_merge_bams.py"
    # shell:
    #     "sambamba merge -t {threads} {output} {input}"


rule mark_duplicates:
    input:
        bam="data/alignments/{sample}.bam",
        ref=get_ref,
    output:
        bam="data/alignments/{sample}_markdup.bam",
    threads: 4
    envmodules:
        config['envmodules']['sambamba']
    shell:
        "sambamba markdup -t {threads} {input.bam} {output.bam}"


rule index_bams:
    input:
        "data/alignments/{sample}_markdup.bam",
    output:
        "data/alignments/{sample}_markdup.bam.bai",
    envmodules:
        config['envmodules']['sambamba']
    shell:
        "sambamba index {input}"
