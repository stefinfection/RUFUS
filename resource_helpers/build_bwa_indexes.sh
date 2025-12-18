#!/bin/bash
host_ref_file="$1"
ref_base=$(basename ${host_ref_file})

# ENV override
: "${RUFUS_ROOT:=/opt/RUFUS}"

bwa="$RUFUS_ROOT/bin/externals/bwa/src/bwa_project/bwa"
BWA_INDEX_CONT_DIR=./rufus_temp/bwa_indexes/ # Needs to be working dir

# Handle .gz extension
if [[ "$ref_base" == *.gz ]]; then
    ref_base_no_gz="${ref_base%.gz}"
else
    ref_base_no_gz="$ref_base"
fi

fasta_idx="$BWA_INDEX_CONT_DIR/${ref_base_no_gz}"

# Index the file
$bwa index -a bwtsw "$fasta_idx"
samtools faidx "$fasta_idx"

echo "These indexes were created by RUFUS for an intermediate BWA step." > $BWA_INDEX_CONT_DIR/README.md
echo "If you wish to reuse them for the next run to save some time, copy them into the same directory as your REFERENCE_FASTA." >> $BWA_INDEX_CONT_DIR/README.md