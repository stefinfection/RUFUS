#!/bin/bash

region=$1
bam=$2
ref=$3

bcftools mpileup -Ov -d 100 -f $ref -r $region $bam > mpileup_${region}.vcf
