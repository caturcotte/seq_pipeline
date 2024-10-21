# ruleorder: merge_bams > mv_nolane_bams

# ignore - moved them all because glob wassn't working
rule symlink_pod5s:
    input:
        find_pod5,
    output:
        "data/pod5/{pod5}.pod5",
    shell:
        "ln -s {input} {output}"


rule make_dorado_sample_sheet:
    input:
        "sample_sheet.csv",
    output:
        "data/resources/dorado_sheet.csv",
    params:
        experiment_id=config["experiment_id"],
    run:
        sheet = pd.read_csv(input[0])
        sheet[[
                "start_time",
                "device_id",
                "position_id",
                "flow_cell_id",
                "protocol_run_id",
            ]]= sheet["ont_folder"].str.split("_", expand=True)
        sheet.to_csv('test2.csv')
        sheet['sample_num'] = sheet['sample_num'].astype(str)
        sheet['sample_num'] = sheet['sample_num'].str.zfill(3)
        sheet['alias'] = sheet[['genotype', 'sample_num']].agg('_'.join, axis=1)
        sheet['barcode'] = sheet['barcode'].astype(str)
        sheet["barcode"] = "barcode" + sheet['barcode'].str.zfill(2)
        sheet["experiment_id"] = params.experiment_id
        sheet = sheet[
            [
                "experiment_id",
                "kit",
                "flow_cell_id",
                "position_id",
                "protocol_run_id",
                "alias",
                "barcode",
            ]
        ]
        sheet.to_csv(output[0], index=False)


# doing alignment at the same time at basecalling takes the same amount of time as basecalling alone
rule ont_basecall_align_demux:
    input:
        dorado=f"workflow/tools/dorado/dorado-0.8.1-{os}-{arch}/bin/dorado",
        #TODO generalize this
        pod5s=[str(i) for i in Path("/users/c/a/cannecar/ndj_progeny/all_pod5s/").glob("*.pod5")],
        ref="data/resources/genome.mmi",
        sample_sheet="data/resources/dorado_sheet.csv",
    output:
        multialn="data/alignments/multialn.bam",
        demux=expand("data/alignments/{progeny}.bam", progeny=all_progeny),
    params:
        in_dir="/users/c/a/cannecar/ndj_progeny/all_pod5s/",
        out_dir=lambda w, output: str(Path(output[1]).parents[0]),
        kit=config["ont_kit"],
        mdl="sup",
    shell:
        "{input.dorado} basecaller {params.mdl} --kit-name \'{params.kit}\' --reference {input.ref} --sample-sheet {input.sample_sheet} {params.in_dir} > {output.multialn} && {input.dorado} demux --kit-name \'{params.kit}\' --sample-sheet={input.sample_sheet} --no-trim --sort-bam --output-dir={params.out_dir} {output.multialn}"


rule align_bowtie2:
    input:
        reads=["data/reads/{parent}_1.fq.gz", "data/reads/{parent}_2.fq.gz"],
        ref=multiext(
            "data/resources/genome",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    output:
        temp("data/alignments/{parent}_raw.bam"),
    params:
        ref_basename="data/resources/genome",
    threads: 32
    envmodules:
        config["envmodules"]["bowtie2"],
    # conda:
    #     "envs/bowtie2.yaml"
    shell:
        "bowtie2 -x {params.ref_basename} -1 {input.reads[0]} -2 {input.reads[1]} -p {threads} | samtools view -1 -o {output}"


rule fix_mate_pairs:
    input:
        "data/alignments/{parent}_raw.bam",
    output:
        temp("data/alignments/{parent}_fixed.bam"),
    envmodules:
        config["envmodules"]["samtools"],
    shell:
        "samtools fixmate -m -O bam,level=1 {input} {output}"


rule sort_bams:
    input:
        "data/alignments/{parent}_fixed.bam",
    output:
        temp("data/alignments/{parent}_sorted.bam"),
    threads: 8
    envmodules:
        config["envmodules"]["samtools"],
    shell:
        "samtools sort {input} -l 1 -o {output} --threads {threads}"


# marking duplicates not necessary for non-PCR based Nanopore kits
rule mark_duplicates:
    input:
        bam="data/alignments/{parent}_sorted.bam",
        ref="data/resources/genome",
    output:
        "data/alignments/{parent}.bam",
    threads: 4
    envmodules:
        config["envmodules"]["sambamba"],
    shell:
        "sambamba markdup -t {threads} {input.bam} {output}"
