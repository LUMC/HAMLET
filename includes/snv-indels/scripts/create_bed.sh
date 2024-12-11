#!/usr/bin/env bash

readonly fasta=$1
readonly gtf=$2

# Get the chromosomes from the fasta file
grep ">" ${fasta} | cut -f 1 -d ' ' | tr -d ">" > chroms.txt

# Create a BED file out of the GTF, for the chromosomes we found before
cat ${gtf} | awk 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,$10,$14,$7}}' | tr -d ';' | grep -f chroms.txt > all_genes.bed

# Merge all regions that overlap, ignoring the strand. For the human genome, this gives 33k regions instead of 61k
sort -k1,1 -k2,2n all_genes.bed | bedtools merge
