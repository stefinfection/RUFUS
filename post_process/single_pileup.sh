#!/bin/bash

# Takes in a list of regions to perform a pileup on in bam, returns a single bgzipped pileup vcf
regions=$1
bam=$2
ref=$3

DEPTH=500
BCFTOOLS="/opt/bcftools/bcftools"

$BCFTOOLS mpileup -Oz -d $DEPTH -f "$ref" -r "${regions}" "$bam"