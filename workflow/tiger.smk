rule convert_tsv_to_parquet:
    input:
        "data/tsvs/{sample}.tsv",
    output:
        "data/parquets/{sample}_converted.parquet",
    conda:
        "envs/convert_tsv_to_parquet.yaml"
    script:
        "scripts/convert_tsv_to_parquet.py"


rule rearrange_reference_parquet:
    input:
        ref=f"data/parquets/{config['ref_parent']}_converted.parquet",
        alt=f"data/parquets/{config['alt_parent']}_converted.parquet",
    output:
        "data/parquets/reference.parquet",
        db=temp("reference.db"),
    params:
        ref_name=config["ref_parent"],
        alt_name=config["alt_parent"],
        min_qual=config["min_qual"],
        min_depth=config["min_depth"],
    # threads: 16
    shell:
        "ipython workflow/scripts/rearrange_reference_parquet.py -- {input.ref} {input.alt} {params.ref_name} {params.alt_name} {params.min_qual} {params.min_depth} {output[0]} {output.db}"


rule rearrange_sample_parquet:
    input:
        smp="data/parquets/{sample}_converted.parquet",
        ref="data/parquets/reference.parquet",
    output:
        "data/parquets/{sample}.parquet",
        temp("{sample}.db"),
    # threads: 16
    shell:
        "ipython workflow/scripts/rearrange_sample_parquet.py -- {wildcards.sample} {input.ref} {input.smp} {output[0]}"


rule make_tiger_input:
    input:
        "data/parquets/{sample}.parquet",
    output:
        "data/tiger/{sample}/{sample}.tiger_input.txt",
    params:
        chr_names=list(config["chromosomes"].keys()),
    run:
        tiger_chrom_name_dict = dict(
            zip(params["chr_names"], range(1, len(params["chr_names"])+1))
        )
        df = pd.read_parquet(input[0])
        df = df.loc[df['chromosome'].isin(params['chr_names'])]
        df = df.replace({"chromosome": tiger_chrom_name_dict})
        df = df.rename(columns={"chromosome": "tiger_chrom"})
        df = df.fillna(0)
        df = df.astype({
            'tiger_chrom': 'str',
            'position': 'int64',
            'ref_reads': 'int32',
            'variant_reads': 'int32'
        })
        df = df[['tiger_chrom', 'position', 'reference', 'ref_reads', 'variant', 'variant_reads']]
        df.to_csv(
            output[0],
            sep="\t",
            header=False,
            index=False,
            columns=[
                "tiger_chrom",
                "position",
                "reference",
                "ref_reads",
                "variant",
                "variant_reads",
            ],
        )


rule make_chrom_lengths_file:
    output:
        "data/resources/chrom_lengths_tiger.txt",
    run:
        chr_nums = range(1, len(list(config['chromosomes'].keys()))+1)
        print(chr_nums)
        with open(output[0], 'w+') as file:
            for key, val in config['chromosomes']:
                file.write(f"{key}\t{val}")


rule run_rtiger:
    input:
        "data/tiger/{sample}/{sample}.tiger_input.txt",
    output:
        "data/tiger/{sample}/GenotypePlot.pdf"
    params:
        sample_name="{sample}",
        generation=lambda w: get_fields(w.sample)['sample_type']
    # conda:
    #     "envs/rtiger.yaml"
    script:
        "scripts/run_rtiger.R"


rule run_basecaller:
    input:
        "data/tiger/{sample}/{sample}.tiger_input.txt",
    output:
        "data/tiger/{sample}/{sample}.base_caller.txt",
    # envmodules:
    #     config["envmodules"]["java"],
    conda:
        "envs/java.yaml"
    shell:
        "java -jar workflow/scripts/TIGER_scripts/base_caller.jar -r {input} -o {output} -n bi"


rule allele_frequency_estimator:
    input:
        "data/tiger/{sample}/{sample}.tiger_input.txt",
    output:
        "data/tiger/{sample}/{sample}.bmm_allele_freqs.txt",
    envmodules:
        config["envmodules"]["java"],
    shell:
        "java -jar workflow/scripts/TIGER_scripts/allele_freq_estimator.jar -r {input} -o {output} -n bi -w 1000"


rule apply_beta_mixture_model:
    input:
        "data/tiger/{sample}/{sample}.bmm_allele_freqs.txt",
    output:
        "data/tiger/{sample}/{sample}.bmm_intersections.txt",
    envmodules:
        config["envmodules"]["r"],
    shell:
        "Rscript --vanilla workflow/scripts/TIGER_scripts/beta_mixture_model.R {input} {output}"


rule prepare_for_hmm_estimation:
    input:
        tiger="data/tiger/{sample}/{sample}.tiger_input.txt",
        basecaller="data/tiger/{sample}/{sample}.base_caller.txt",
        chrom_lengths="data/resources/dm6_chrom_lengths_tiger.txt",
    output:
        "data/tiger/{sample}/{sample}.tiger.hmm_prep_prob.txt",
    envmodules:
        config["envmodules"]["perl"],
    shell:
        "perl workflow/scripts/TIGER_scripts/prep_hmm_probs.pl -s {wildcards.sample} -m {input.tiger} -b {input.basecaller} -c {input.chrom_lengths} -o {output}"


rule calculate_transmission_emission:
    input:
        allele_freqs="data/tiger/{sample}/{sample}.bmm_allele_freqs.txt",
        prep="data/tiger/{sample}/{sample}.tiger.hmm_prep_prob.txt",
        isecs="data/tiger/{sample}/{sample}.bmm_intersections.txt",
        chrom_lengths="data/resources/dm6_chrom_lengths_tiger.txt",
    output:
        "data/tiger/{sample}/{sample}_hmm_model",
    envmodules:
        config["envmodules"]["perl"],
    shell:
        "perl workflow/scripts/TIGER_scripts/hmm_probs.pl -s {input.allele_freqs} -p {input.prep} -a {input.isecs} -c {input.chrom_lengths} -o data/tiger/{wildcards.sample}/{wildcards.sample}"


rule run_hmm:
    input:
        basecaller="data/tiger/{sample}/{sample}.base_caller.txt",
        model="data/tiger/{sample}/{sample}_hmm_model",
    output:
        "data/tiger/{sample}/{sample}.hmm_output.txt",
    envmodules:
        config["envmodules"]["java"],
    shell:
        "java -jar workflow/scripts/TIGER_scripts/hmm_play.jar -r {input.basecaller} -z {input.model} -t bi -o {output}"


rule estimate_breakpoints:
    input:
        tiger="data/tiger/{sample}/{sample}.tiger_input.txt",
        hmm="data/tiger/{sample}/{sample}.hmm_output.txt",
        chrom_lengths="data/resources/dm6_chrom_lengths_tiger.txt",
    output:
        "data/tiger/{sample}/{sample}.CO_estimates.txt",
    envmodules:
        config["envmodules"]["perl"],
    shell:
        "perl workflow/scripts/TIGER_scripts/estimate_breaks.pl -s {wildcards.sample} -m {input.tiger} -b {input.hmm} -c {input.chrom_lengths} -o {output}"


rule refine_breakpoints:
    input:
        tiger="data/tiger/{sample}/{sample}.tiger_input.txt",
        cos="data/tiger/{sample}/{sample}.CO_estimates.txt",
    output:
        "data/tiger/{sample}/{sample}.CO_estimates.breaks.txt",
    envmodules:
        config["envmodules"]["perl"],
    shell:
        "touch {output} && perl workflow/scripts/TIGER_scripts/refine_recombination_break.pl {input.tiger} {output}"


rule smooth_breakpoints:
    input:
        "data/tiger/{sample}/{sample}.CO_estimates.breaks.txt",
    output:
        "data/tiger/{sample}/{sample}.smoothed.breaks.txt",
    envmodules:
        config["envmodules"]["perl"],
    shell:
        "perl workflow/scripts/TIGER_scripts/breaks_smoother.pl -b {input} -o {output}"


rule process_tiger_output:
    input:
        co_estimates=aggregate_input,
        breaks=lambda w: aggregate_input(w, breaks=True),
        # co_estimates=expand(
        #     "data/tiger/{sample}/{sample}.CO_estimates.txt",
        #     sample=progeny),
        # breaks=expand(
        #     "data/tiger/{sample}/{sample}.CO_estimates.breaks.txt", sample=progeny
        #     ),
    output:
        hmm_states="data/tiger_output/hmm_states.csv",
        hmm_intervals="data/tiger_output/hmm_intervals.csv"
    conda:
        "envs/process_tiger_outputs.yaml"
    threads: 64
    script:
        "scripts/tiger_output.py"
