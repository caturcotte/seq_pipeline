#!/bin/bash

module load samtools

bcftools view -v ref,snps "${1}" | bcftools query -f \
	'[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%GT\t%DP\t%AD]\n' \
	-o "${2}" &&
	sed -i '1s/^/sample\tchromosome\tposition\treference\tvariant\tquality\tgenotype\tdepth\tallele_depth\n/' "${2}" 
