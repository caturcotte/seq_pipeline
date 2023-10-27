rule bwa_mem:
    input:
        r1="reads/{sample}_1.fq.gz",
        r2="reads/{sample}_2.fq.gz",
        ref=get_ref,
        ref_idx=lambda w: get_ref(w, fai=True),
    output:
        sam="mapped/{sample}_bwa.sam",
    resources:
        time="5-0",
        mem_mb=15000,
    threads: 8
    shell:
        "bwa mem {input.ref} {input.r1} {input.r2} > {output}"


rule bowtie2:
    input:
        r1="reads/{sample}_1.fq.gz",
        r2="reads/{sample}_2.fq.gz",
        ref=get_ref_bowtie2,
    output:
        sam="mapped/{sample}_bt2.sam",
    params:
        extra=lambda w: f"--rg-id={w.sample} --rg SM:{w.sample}",
        ref_basename=get_ref,
    resources:
        time="5-0",
        mem_mb=15000,
    threads: 8
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
        mem_mb=15000,
    threads: 8
    shell:
        "minimap2 -ax map-ont {input.ref} {input.reads} --threads {threads} > {output}"


rule fix_mate_pairs:
    input:
        get_aligned_reads,
    output:
        temp("mapped/{sample}.bam"),
    resources:
        time="2:00:00",
        mem_mb=15000,
    shell:
        "samtools fixmate -m -z on -O bam,level=1 {input} {output}"


rule samtools_sort:
    input:
        "mapped/{sample}.bam",
    output:
        temp("mapped/{sample}_sort.bam"),
    threads: 8
    resources:
        time="2:00:00",
        mem_mb=15000,
    shell:
        "samtools sort {input} -l 1 -o {output} --threads {threads}"


rule mark_duplicates:
    input:
        bam="mapped/{sample}_sort.bam",
        ref=get_ref,
    output:
        bam="mapped/{sample}_sort_dedup.bam",
        metrics="metrics/{sample}_dedup_metrics.txt",
    threads: 8
    resources:
        time="1-0",
        mem_mb=50000,
    params:
        d=get_opt_dup_distance,
    shell:
        "samtools markdup -d {params.d} -r -f {output.metrics} --threads {threads} {input.bam} {output.bam}"


rule samtools_index:
    input:
        "mapped/{sample}_sort_dedup.bam",
    output:
        "mapped/{sample}_sort_dedup.bam.bai",
    threads: 4
    resources:
        time="1-0",
    wrapper:
        "v2.3.1/bio/samtools/index"
