#!/bin/bash

host_ref_file="$1"
ref_base=$(basename ${host_ref_file})
fasta_idx="/mnt/${ref_base}"
$bwa index -a bwtsw $fasta_idx
$samtools faidx $fasta_idx