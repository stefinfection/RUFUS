#!/bin/bash

chr=$1
start_coord=$2
end_coord=$3
bam=$4
ref=$5

echo "$bam and $ref" >&2
bcftools mpileup -Ov -d 100 -f $ref -r "${chr}:${start_coord}-${end_coord}" $bam > mpileup_${region}.vcf
