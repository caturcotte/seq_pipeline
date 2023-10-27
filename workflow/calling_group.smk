rule separate_into_samples:
    input:
        get_group_from_sample,
    output:
        "called/{sample}_{caller}_unprocessed.bcf"
    resources:
        time="2:00:00",
    shell:
        "bcftools view -s {wildcards.sample} -a -o {output} {input}"

