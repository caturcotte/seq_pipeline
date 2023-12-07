rule align_bwa:
    input:
        reads=get_reads_to_map,
        ref=get_ref,
        ref_bwa_idx=get_ref_bwa,
    output:
        "data/alignments/{sample}_{iden}_bwa.bam",
    resources:
        time="5-0",
    threads: 16
    conda:
        "envs/bwa.yaml"
    shell:
        "bwa mem {input.ref} {input.reads} -R '@RG\\tID:{wildcards.iden}\\tSM:{wildcards.sample}' --threads {threads} | samtools view -1 -o {output}"


rule align_bowtie2:
    input:
        reads=get_reads_to_map,
        ref=get_ref_bowtie2,
    output:
        "data/alignments/{sample}_{iden}_bt2.bam",
    params:
        ref_basename=get_ref,
    resources:
        time="5-0",
    threads: 16
    conda:
        "envs/bowtie2.yaml"
    shell:
        "bowtie2 -x {params.ref_basename} -1 {input.reads[0]} -2 {input.reads[1]} -p {threads} --rg-id {wildcards.iden} --rg 'SM:{wildcards.sample}' | samtools view -1 -o {output}"


rule align_minimap2:
    input:
        reads=get_reads_to_map,
        ref=get_ref_minimap2,
    output:
        "data/alignments/{sample}_{iden}_mm2.bam",
    resources:
        time="5-0",
    threads: 16
    conda:
        "envs/minimap2.yaml"
    shell:
        "minimap2 -ax map-ont {input.ref} {input.reads} -R '@RG\\tID:{wildcards.iden}\\tSM:{wildcards.sample}' --threads {threads} | samtools view -1 -o {output}"


rule fix_mate_pairs:
    input:
        get_aligned_reads,
    output:
        temp("data/alignments/{sample}_{iden}.bam"),
    resources:
        time="2:00:00",
    shell:
        "samtools fixmate -m -O bam,level=1 {input} {output}"


rule sort_bams:
    input:
        "data/alignments/{sample}_{iden}.bam",
    output:
        temp("data/alignments/{sample}_{iden}_sort.bam"),
    threads: 8
    resources:
        time="2:00:00",
        mem_mb=20000,
    shell:
        "samtools sort {input} -l 1 -o {output} --threads {threads}"


rule merge_bams:
    input:
        get_alns_to_merge,
    output:
        "data/alignments/{sample}.bam",
    threads: 4
    shell:
        "sambamba merge -t {threads} {output} {input}"


rule remove_duplicates:
    input:
        bam="data/alignments/{sample}.bam",
        ref=get_ref,
    output:
        bam="data/alignments/{sample}_dedup.bam",
    threads: 4
    resources:
        time="1-0",
    log:
        "data/qc/sambamba/{sample}.log",
    shell:
        "sambamba markdup -r -t {threads} {input.bam} {output.bam} > {log} 2>&1"


rule index_bams:
    input:
        "data/alignments/{sample}_dedup.bam",
    output:
        "data/alignments/{sample}_dedup.bam.bai",
    resources:
        time="2:00:00",
    shell:
        "sambamba index {input}"
