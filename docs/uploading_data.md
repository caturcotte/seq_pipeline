---
layout: page
title: "Uploading data to Longleaf"
permalink: /uploading-data
---

## Uploading Novogene data

When your data is ready, you should get an email that looks like this. Click the link highlighted in yellow:

![[Pasted image 20231120164637.png]]

It should ask you for your email and password (you used these to make an account to put in the order).  

Go to the data release for your project and click Batch download. A panel will open up on the side with a [wget](https://www.gnu.org/software/wget/) command. In the terminal, connect to Longleaf and go to the directory in which you plan to store your data. Send the wget command as a job to [SLURM](https://help.rc.unc.edu/longleaf-slurm-examples/), Longleaf's job handling application. This will take a while.

```
# connect to Longleaf
$ ssh <onyen>@longleaf.unc.edu

# go to the directory where you want to store the data
$ cd /proj/sekellab/projects/<PROJECT-NAME>/

# submit the wget command as a job to SLURM
$ sbatch -p general -N 1 -n 1 --mem=5g -t 05-00:00:00 --wrap="<WGET-COMMAND-FROM-NOVOGENE>"
```


You should be in the directory with your Novogene data, and your directory should now look something like this (with sample being a placeholder for a sample name):

```
$ tree
.
└── usftp21.novogene.com/
    ├── 01.RawData/
    │   └── sample/
    │       ├── MD5.txt
    │       ├── sample_CKDN230029067-1A_HCTF3DSX7_L3_1.fq.gz
    │       ├── sample_CKDN230029067-1A_HCTF3DSX7_L3_2.fq.gz
    │       ├── sample_CKDN230029067-1A_HGKT2DSX7_L4_1.fq.gz
    │       └── sample_CKDN230029067-1A_HGKT2DSX7_L4_2.fq.gz
    ├── 02.Report_X202SC23082820-Z01-F001.zip
    ├── checkSize.xls
    ├── MD5.txt
    ├── MD5.zip
    └── Readme.html
```

>[!note]
The number of files and the suffix for each of the sample files (e.g. CKDN230029067-1A_HCTF3DSX7) will vary.

> [!tip]
> You can use the command `tree` to display a directory structure like the one above.

### Verifying that your download worked

You should check the [[checksum]] of all of your files to make sure the downloads were complete.
```
$ cd usftp21.novogene.com
$ md5sum -c MD5.txt
```

If it works, it will say something like:
```
01.RawData/sample/sample_CKDN230029067-1A_HCTF3DSX7_L3_1.fq.gz: OK
# ... and so on for every file in all of the sample directories
```

>[!failure]
If the md5sum check fails, your data is corrupted and you need to redownload it.

## Uploading Nanopore data

Because the file size to transfer is large, you will need to use [GLOBUS](https://help.rc.unc.edu/getting-started-with-globus-connect) to transfer the data from the Nanopore computer to Longleaf.

First, compress the data to make it easier to transfer. Because /var isn't a subdirectory of your home directory, you will need to perform this operation as the root user (or "super user", `sudo` = **S**uper **U**ser **Do**. This is like running a program with admin privileges, in Windows). Since you are performing this command as root, it will prompt you for the computer's password.

```
$ cd /var/lib/minknow/data/<PROJECT-NAME>/<FLOWCELL-ID>
$ sudo tar -czvf <PROJECT-NAME>.tar.gz fastq_pass/*
```

You should now have a file called `PROJECT-NAME.tar.gz` (where PROJECT-NAME is a placeholder for whatever you named your project in MinKNOW).

GLOBUS can't recognize any directories above your home directory, so let's move the compressed file somewhere else:
```
$ sudo mv <PROJECT-NAME>.tar.gz ~
```

Now `PROJECT-NAME.tar.gz` is in our home directory (`~`).

Start up Globus Connect Personal on the Nanopore computer:
```
$ ~/globusconnectpersonal & disown
```

A window should pop up. Click "Connect" and make sure not to close the window.

In a web browser, go to https://www.globus.org and click "Log in" at the top right. Use the organizational login, and search for University of North Carolina - Chapel Hill in the drop down. There will be a few screens to go through and then you'll be prompted for your ONYEN. See https://help.rc.unc.edu/getting-started-with-globus-connect for more detail and screenshots.

Once you're logged in, click the File Manager tab on the left side. You'll see a split screen view of two different data collections:

![[nanopore.png]]

In the left collection search box, type "nanopc" and click the result that pops up. This is the Nanopore computer. In the right box, search "UNC, Research Computing, DataMover" and click the result that pops up. This is Longleaf.

>[!failure]
>If the icon for nanopc is red, indicating that the connection is offline, start Globus Connect Personal again by typing `~/globusconnectpersonal` in the terminal.

You will automatically be in the home directory for both collections. In the UNC DataMover collection, change the path to `/work/users/<o>/<n>/<onyen>/`.

Find and click your `<PROJECT-NAME>.tar.gz` file in the left collection (it should be highlighted in blue), then click the blue Start button on the left. A green popup should appear saying the transfer request has been submitted.
![[transfer.png]]

Now click the Activity panel on the left side menu. You should find an entry that says "nanopc to UNC, Research Computing, DataMover". Click it,and it will take you to a page showing the progress of your transfer.

When the condition says "SUCCEEDED" your file transfer has finished.

## Cloning the pipeline

Go to your work directory on Longleaf at `/work/users/o/n/onyen/` (where 'o' and 'n' are the first and second letters of your ONYEN, respectively).

```
$ cd /work/users/<o>/<n>/<onyen>/
```

Clone the pipeline from [here](https://github.com/caturcotte/seq_pipeline/tree/main) (click the green button that says Code, then copy the link in the drop down), then go to the directory of the repository.
```
$ git clone https://github.com/caturcotte/seq_pipeline.git
$ cd seq_pipeline/
```

The directory tree of the pipeline should look like this:
```
$ tree
.
├── config.yaml
├── miniforge_installer.sh
├── README.md
├── sample_sheet.csv
├── slurm
│   └── config.yaml
├── slurm_sub.sh
├── Snakefile
└── workflow
    ├── alignment.smk
    ├── calling.smk
    ├── envs
    │   ├── bcftools.yaml
    │   ├── bowtie2.yaml
    │   ├── bwa.yaml
    │   ├── concat_vcfs.yaml
    │   ├── environment.yaml
    │   ├── freebayes.yaml
    │   ├── minimap2.yaml
    │   ├── repeatmasker.yaml
    │   └── samtools.yaml
    ├── functions.smk
    ├── misc.smk
    ├── qc.smk
    └── vcf_filtering.smk
```

We are now ready to [[Running the DNA sequencing analysis pipeline#Cloning the pipeline|set up the pipeline to run]].
