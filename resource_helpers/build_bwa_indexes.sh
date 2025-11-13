#!/bin/bash
host_ref_file="$1"
ref_base=$(basename ${host_ref_file})

bwa="/opt/RUFUS/bin/externals/bwa/src/bwa_project/bwa"
samtools=/opt/samtools/samtools

# Handle .gz extension
if [[ "$ref_base" == *.gz ]]; then
    ref_base_no_gz="${ref_base%.gz}"
else
    ref_base_no_gz="$ref_base"
fi

fasta_idx="/mnt/bwa_indexes/${ref_base}"

# Create bwa_indexes directory if it doesn't exist
mkdir -p /mnt/bwa_indexes

# Index the file
$bwa index -a bwtsw "$fasta_idx"
$samtools faidx "$fasta_idx"

echo "These indexes were created by RUFUS for an intermediate BWA step." > /mnt/bwa_indexes/README.md
echo "If you wish to reuse them for the next run to save some time, copy them into the same directory as your REFERENCE_FASTA." >> /mnt/bwa_indexes/README.md