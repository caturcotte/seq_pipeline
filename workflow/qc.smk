# general QC checks for reads
rule fastqc_single:
    input:
        "data/reads/{sample}.fq.gz",
    output:
        html="data/qc/fastqc/{sample}.html",
        zip="data/qc/fastqc/{sample}_fastqc.zip",
    wrapper:
        "v2.11.1/bio/fastqc"


rule fastqc_paired:
    input:
        "data/reads/{sample}_{read}.fq.gz",
    output:
        html="data/qc/fastqc/{sample}_{read}.html",
        zip="data/qc/fastqc/{sample}_{read}_fastqc.zip",
    wrapper:
        "v2.11.1/bio/fastqc"


# read depth calculations
rule mosdepth:
    input:
        bam="data/alignments/{sample}_dedup.bam",
        bai="data/alignments/{sample}_dedup.bam.bai",
    output:
        "data/qc/mosdepth/{sample}.mosdepth.global.dist.txt",
        "data/qc/mosdepth/{sample}.per-base.bed.gz",
        summary="data/qc/mosdepth/{sample}.mosdepth.summary.txt",
    threads: 4
    wrapper:
        "v2.11.1/bio/mosdepth"


# rule fastq_screen_get_genomes:
#   output:
#     "data/resources/FastQ_Screen_Genomes/fastq_screen.conf"
#   conda:
#     "envs/fastq_screen.yaml"
#   shell:
#     "fastq_screen --get_genomes --outdir data/resources/"

# # detect contamination from other species
# rule fastq_screen:
#   input:
#     "data/reads/{sample}.fq",
#   output:
#     txt="data/qc/{sample}.fastq_screen.txt",
#     png="data/qc/{sample}.fastq_screen.png",
#   params:
#     fastq_screen_config="data/resources/FastQ_Screen_Genomes/fastq_screen.conf"
#     subset=100000,
#     aligner="bwa",
#   threads: 8
#   wrapper:
#     "v2.11.1/bio/fastq_screen"
# rule gatk_create_dict:
#   input:
#     "resources/{ref}.fa",
#   output:
#     "resources/{ref}.dict"
#   wrapper:
#     "v2.13.0/bio/picard/createsequencedictionary"

# # stats for vcf files
# rule gatk_varianteval:
#   input:
#     vcf="called/{sample}_{caller}_norm.bcf",
#     # vcf="data/calls/{sample}.vcf",
#     ref=get_ref,
#     dict="genome.dict",
#   output:
#     vcf="{sample}.varianteval.grp",
#   wrapper:
#     "v2.11.1/bio/gatk/varianteval"

# detects potential sample mixture, and sanity checks whether samples are derived from closely related individuals (sibling/parents)
# rule verify_bam_id:
#   input:
#     bam="data/alignments/{sample}_processed.bam",
#     ref=get_ref,
#   output:
#     selfsm="data/qc/verify_bam_id/{sample}.selfSM",
#     ancestry="data/qc/verify_bam_id/{sample}.ancestry",
#   wrapper:
#     "v2.11.1/bio/verifybamid/verifybamid2"


rule vcf_stats:
    input:
        calls=get_caller,
        # ref=get_ref,
    output:
        "data/qc/bcftools/{sample}_{caller}.stats.txt",
    # conda:
    #     "envs/bcftools.yaml"
    wrapper:
        "v2.13.0/bio/bcftools/stats"


rule multiqc:
    input:
        # expand("data/qc/verify_bam_id/{sample}.selfSM", sample=samples),
        expand(
            "data/qc/bcftools/{sample}_{caller}.stats.txt",
            sample=samples,
            caller=config["caller"],
        ),
        glob_wildcards("data/qc/cutadapt/{sample}_{prefix}.txt"),
        expand("data/qc/sambamba/{sample}.log", sample=samples),
        # expand("data/qc/varianteval/{sample}.varianteval.grp", sample=samples),
        # expand("data/qc/fastq_screen/{sample}.fastq_screen.txt", sample=samples),
        expand("data/qc/fastqc/{sample}_{read}_fastqc.zip", sample=samples, read=[1, 2]),
        expand("data/qc/mosdepth/{sample}.mosdepth.global.dist.txt", sample=samples),
    output:
        "data/qc/multiqc.html",
        directory("data/qc/multiqc_data"),
    params:
        extra="--data-dir",
    wrapper:
        "v2.11.1/bio/multiqc"
