---
layout: page
title: "Understanding the pipeline"
permalink: /understanding-the-pipeline
---

## How the pipeline works

#### Alignment and [BAM](https://en.wikipedia.org/wiki/Binary_Alignment_Map) processing

First, all of the reads or read pairs for a given sample will be individually aligned to whatever reference genome is selected in `config.yaml`. The alignments will then be sorted and merged into a single BAM file, with lane information being preserved in the [read group](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups) tags of the reads. Next, the BAM files will be processed to remove [duplicates](https://chaochungkuo.github.io/notes/optical-nonoptical-duplicates/). This is only necessary for Illumina sequencing and sequencing of PCR-amplified DNA, but it will not be detrimental for analysis of Nanopore reads.
```mermaid
flowchart TD;
	A[sample 1]--- a{{sequencing}}-->B[sample 1 lane 1 reads] & C[sample 1 lane 2 reads] & D[sample 1 lane 3 reads]
	B --- b{{alignment}} -->E[sample 1 lane 1 BAM]
	C --- c{{alignment}} -->F[sample 1 lane 2 BAM]
	D --- d{{alignment}} -->G[sample 1 lane 3 BAM]
	E & F & G --> q{{sort and merge}} --> H[sample 1 BAM]
	H --- f{{remove duplicates}} --> I[sample 1 BAM, duplicates removed]
```

#### Group variant calling

Next, the polished BAM files (with duplicates removed) will be pooled into the groups specified in the sample sheet. Variants will be called and sorted into [[Variant Call Format]] files.
```mermaid
flowchart TD
	H[sample 1 BAM] & I[sample 2 BAM] & J[sample 3 BAM]---z{group 1}---e{{variant calling}} -->K[VCF]
	K --- a{{filtering}} --> L[filtered VCF] ---r{{separate by sample}} -->M[sample 1 VCF] & N[sample 2 VCF] & O[sample 3 VCF]
	P[sample 4 BAM] & Q[sample 5 BAM]---g{group 2}-->k[the same process, separately]
```

Using these VCFs, you can make consensus sequences for your samples by applying the variants in the sample to the reference genome. The resulting FASTA files can be viewed in Snapgene, etc. for downstream analysis.
```mermaid
flowchart TD
	A[sample 1 VCF] & B[reference FASTA] ---q{{make consensus sequence}}-->C[sample 1 reference FASTA]
```

