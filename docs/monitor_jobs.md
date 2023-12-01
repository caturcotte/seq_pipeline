---
layout: page
title: "Monitor the status of your jobs"
permalink: /monitor-jobs
---

# Monitor the status of all of your jobs

You can do the following to monitor the status of your jobs:
```
$ squeue -u <onyen>
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
```

The pipeline will submit jobs to SLURM automatically through your account, so you will likely see several jobs here. ST indicates the status of the job, where PD means pending and R means that the job is currently running. TIME indicates the time that the job has been running for.

If you want to continually monitor the status of your jobs:
```
$ squeue -u <onyen> -i 60
```

This will update with the status of your jobs every minute (press Ctrl+C to cancel).

## Monitor the output log

To continually follow the output of your jobs, use:
```
$ tail -F slurm-<JOB-NUMBER>.out
```

This will continually report the status of all of the jobs (press Ctrl+C to cancel).

If you see an error, such as:
```
[Mon Nov 27 10:41:49 2023]
Error in rule sort_bams:
    jobid: 15
    input: data/alignments/ndj_01_CKDN230029052-1A_HGL2TDSX7_L1.bam
    output: data/alignments/ndj_01_CKDN230029052-1A_HGL2TDSX7_L1_sort.bam
    shell:
        samtools sort data/alignments/ndj_01_CKDN230029052-1A_HGL2TDSX7_L1.bam -l 1 -o data/alignments/ndj_01_CKDN230029052-1A_HGL2TDSX7_L1_sort.bam --threads 4
        (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)
    cluster_jobid: Submitted batch job 23700119

Error executing rule sort_bams on cluster (jobid: 15, external: Submitted batch job 23700119, jobscript: /work/users/c/a/cannecar/seq_pipeline/.snakemake/tmp.x1vvuyv2/snakejob.sort_bams.15.sh). For error details see the cluster log and the log files of the involved rule(s).
```

Then there was a problem with a job. You can view the output of one of the jobs that failed to see more specific details (see below).
## View the output of specific jobs

>[!tip]
>If one of the jobs gave you an error. you can use the error message from the main output log to quickly find the corresponding log file.
>
Let's say the error message from the main output logis "Error executing rule format_bcf on cluster (jobid: 26, external: Submitted batch job **23966205**".
>
In seq_pipeline/ you can find the corresponding log file with `find logs/ -name '23966205'` (which will print the name of the file to the terminal output).

In the logs/ folder, there will be folders for each job that has been run. For example:
```
$ cd /work/users/c/a/cannecar/seq_pipeline/logs/
$ ls
align_bwa		bwa_idx	    fix_mate_pairs   mask_repeats  separate_into_samples
bcftools_call		bwa_mem	    format_bcf       merge_bams    sort_bams
bcftools_index		faidx_ref	    freebayes	     merge_tsvs    vcf_stats
bcftools_mpileup_single  fastqc_paired	    index_bams       mosdepth
bcftools_norm		filter_low_quality  mark_duplicates  multiqc
```

If we go into one of these folders, there will be a file for each sample. named in the form SLURM-JOB-ID-JOB-NAME-sample=SAMPLE.out:
```
$ cd merge_bams
$ ls
'merge_bams-sample=ndj_01.out'	'merge_bams-sample=ndj_07.out'	'merge_bams-sample=ndj_13.out'
'merge_bams-sample=ndj_02.out'	'merge_bams-sample=ndj_08.out'	'merge_bams-sample=ndj_14.out'
'merge_bams-sample=ndj_03.out'	'merge_bams-sample=ndj_09.out'	'merge_bams-sample=ndj_15.out'
'merge_bams-sample=ndj_04.out'	'merge_bams-sample=ndj_10.out'	'merge_bams-sample=ndj_16.out'
'merge_bams-sample=ndj_05.out'	'merge_bams-sample=ndj_11.out'
'merge_bams-sample=ndj_06.out'	'merge_bams-sample=ndj_12.out'
```


You can view the output of one of these files using `cat FILE-NAME` (which prints the contents to stdout, the terminal output) or `less FILE-NAME` (which opens the contents in a separate window - press `q` to return to the terminal prompt. If you're wondering why it's called  less, it's because `less` is the improved successor of a more suitably named program called `more`, which does the same thing... so "less is more").

If the job completed successfully, the end of the file will probably say something like
```
[Tue Nov 28 11:26:26 2023]
Finished job 0.
1 of 1 steps (100%) done
```

If not, the job is either still running or there will be some error you can see in the log file.
