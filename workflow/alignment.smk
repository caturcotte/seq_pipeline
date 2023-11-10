rule bwa_mem:
    input:
        r1=lambda w: get_reads(w, r=1),
        r2=lambda w: get_reads(w, r=2),
        ref=get_ref,
        ref_bwa_idx=get_ref_bwa,
    output:
        sam="mapped/{sample}_bwa.sam",
    resources:
        time="5-0",
    threads: 16
    conda:
        "envs/bwa.yaml"
    shell:
        "bwa mem {input.ref} {input.r1} {input.r2} > {output}"


rule bowtie2:
    input:
        r1=lambda w: get_reads(w, r=1),
        r2=lambda w: get_reads(w, r=2),
        ref=get_ref_bowtie2,
    output:
        sam="mapped/{sample}_bt2.sam",
    params:
        extra=lambda w: f"--rg-id={w.sample} --rg SM:{w.sample}",
        ref_basename=get_ref,
    resources:
        time="5-0",
    threads: 16
    conda:
        "envs/bowtie2.yaml"
    shell:
        "bowtie2 -x {params.ref_basename} -1 {input.r1} -2 {input.r2} -p {threads} > {output.sam}"


rule minimap2:
    input:
        reads="reads/{sample}.fq.gz",
        ref=get_ref_minimap2,
    output:
        "mapped/{sample}_mm2.sam",
    resources:
        time="5-0",
    threads: 16
    conda:
        "envs/minimap2.yaml"
    shell:
        "minimap2 -ax map-ont {input.ref} {input.reads} --threads {threads} > {output}"


rule fix_mate_pairs:
    input:
        get_aligned_reads,
    output:
        temp("mapped/{sample}.bam"),
    resources:
        time="2:00:00",
    # conda:
    #     "envs/samtools.yaml"
    shell:
        "samtools fixmate -m -O bam,level=1 {input} {output}"


rule samtools_sort:
    input:
        "mapped/{sample}.bam",
    output:
        temp("mapped/{sample}_sort.bam"),
    threads: 4
    resources:
        time="2:00:00",
    # conda:
    #     "envs/samtools.yaml"
    shell:
        "samtools sort {input} -l 1 -o {output} --threads {threads}"


rule mark_duplicates:
    input:
        bam="mapped/{sample}_sort.bam",
        ref=get_ref,
    output:
        bam="mapped/{sample}_sort_dedup.bam",
        # metrics="metrics/{sample}_dedup_metrics.txt",
    threads: 4
    resources:
        time="1-0",
    # params:
    #     d=get_opt_dup_distance,
    # conda:
    #     "envs/samtools.yaml"
    shell:
        "sambamba markdup -r -t {threads} {input.bam} {output.bam}"


rule samtools_index:
    input:
        "mapped/{sample}_sort_dedup.bam",
    output:
        "mapped/{sample}_sort_dedup.bam.bai",
    resources:
        time="2:00:00",
    # conda:
    #     "envs/samtools.yaml"
    shell:
        "sambamba index {input}"
