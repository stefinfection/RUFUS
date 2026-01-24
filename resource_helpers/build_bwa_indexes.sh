#!/bin/bash
host_ref_file="$1"

# ENV override
: "${RUFUS_ROOT:=/opt/RUFUS}"

bwa="$RUFUS_ROOT/bin/externals/bwa/src/bwa_project/bwa"

# Handle .gz extension
if [[ "$host_ref_file" == *.gz ]]; then
    host_ref_file_no_gz="${ref_base%.gz}"
else
    host_ref_file_no_gz="$ref_base"
fi

fasta_idx="${host_ref_file_no_gz}"

# Index the file
$bwa index -a bwtsw "$fasta_idx"
samtools faidx "$fasta_idx"

echo "These indexes were created by RUFUS for an intermediate BWA step." > README.md
echo "If you wish to reuse them for the next run to save some time, copy them into the same directory as your REFERENCE_FASTA." >> README.md