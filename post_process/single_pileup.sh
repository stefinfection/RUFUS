#!/bin/bash

# Takes in a list of regions to perform a pileup on in bam, returns a single bgzipped pileup vcf
regions=$1
bam=$2
ref=$3

DEPTH=100

if [ ! -z "$start_coord" ] && [ ! -z "$end_coord" ]; then
  bcftools mpileup -Ov -d $DEPTH -f $ref -r "${chr}:${start_coord}-${end_coord}" -o mpileup_${chr}_${start_coord}_${end_coord}.vcf $bam
else
  bcftools mpileup -Ov -d $DEPTH -f $ref -r "${chr}" -o mpileup_${chr}.vcf $bam
fi