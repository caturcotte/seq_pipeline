# general QC checks for reads
rule fastqc_single:
    input:
        "data/reads/{sample}_{iden}{read}.fq.gz",
    output:
        html="data/qc/fastqc/{sample}_{iden}{read}.html",
        zip="data/qc/fastqc/{sample}_{iden}{read}_fastqc.zip",
    wrapper:
        "v2.11.1/bio/fastqc"


# read depth calculations
rule mosdepth:
    input:
        bam="data/alignments/{sample}_markdup.bam",
        bai="data/alignments/{sample}_markdup.bam.bai",
    output:
        "data/qc/mosdepth/{sample}.mosdepth.global.dist.txt",
        "data/qc/mosdepth/{sample}.per-base.bed.gz",
        summary="data/qc/mosdepth/{sample}.mosdepth.summary.txt",
    threads: 4
    wrapper:
        "v2.11.1/bio/mosdepth"


checkpoint make_mosdepth_df:
    input:
        expand("data/qc/mosdepth/{sample}.mosdepth.summary.txt", sample=progeny),
    output:
        "data/qc/mosdepth/all_progeny.csv"
    run:
        df_list = []
        for i in input:
            df = pd.read_csv(i, sep='\t')
            df['sample'] = i.rstrip(".mosdepth.summary.txt").lstrip("data/qc/mosdepth/")
            df_list.append(df)
        dfs = pd.concat(df_list)
        samples_pass = dfs.query("(chrom == 'chr2L' | chrom == 'chr2R') & mean >= 5")
        samples_unique = list(samples_pass['sample'].unique())
        with open(output[0], 'w+') as file:
            for i in samples_unique:
                file.write(f'{i}\n')


# detects potential sample mixture, and sanity checks whether samples are derived from closely related individuals (sibling/parents)
# rule verify_bam_id:
#   input:
#     ref_vcf="data/vcfs/{config['ref_parent']}.vcf.gz"
#     ref_genome=
#     bam="data/alignments/{sample}_processed.bam",
#     ref=get_ref,
#   output:
#     selfsm="data/qc/verify_bam_id/{sample}.selfSM",
#     ancestry="data/qc/verify_bam_id/{sample}.ancestry",
#   wrapper:
#     "v2.11.1/bio/verifybamid/verifybamid2"


rule vcf_stats:
    input:
        calls="data/calls/{sample}_raw.bcf"
        # ref=get_ref,
    output:
        "data/qc/bcftools/{sample}.stats.txt",
    # conda:
    #     "envs/bcftools.yaml"
    wrapper:
        "v2.13.0/bio/bcftools/stats"


# idens = [j for i in samples for j in get_ids_for_sample(i)]


rule multiqc:
    input:
        # expand("data/qc/verify_bam_id/{sample}.selfSM", sample=samples),
        # expand(
        #     "data/qc/bcftools/{sample}_{caller}.stats.txt",
        #     sample=samples,
        #     caller=config["calling"]["caller"],
        # ),
        # expand("data/qc/sambamba/{sample}.log", sample=samples),
        glob_wildcards("data/qc/fastqc/{sample}_{iden}{read}_fastqc.zip"),
        expand("data/qc/mosdepth/{sample}.mosdepth.global.dist.txt", sample=samples),
    output:
        "data/qc/multiqc.html",
        directory("data/qc/multiqc_data"),
    params:
        extra="--data-dir",
    wrapper:
        "v2.11.1/bio/multiqc"
