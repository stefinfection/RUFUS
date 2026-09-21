#!/bin/bash
# Build a synthetic "RUFUS.interpret output" VCF for the finalize_vcf.sh tests.
#
# REF bases are read from the fixture reference so the records are genuinely valid -- otherwise
# `bcftools +fill-from-fasta -c REF` silently rewrites them mid-stage and the test measures nothing.
#
# The planted set is chosen to exercise the representation logic specifically:
#   20000  AAC>CTG   MNP3 whose FIRST atom (A>C at 20000) is the germline SNV carried by
#                    fixtures/somatic/normal.bam -- the PARTIAL-MATCH case that drives the
#                    all-atoms rule in issue #98. Atoms 2 and 3 are novel.
#   60000  T>A       plain SNV, novel
#   70000  ATG>CAC   MNP3, fully novel -- decompose_blocksub splits this, `norm -a` also splits it,
#                    so it is the record that exposes control-path divergence once atomize is gone.
#   80000  G>GTTTT   insertion, novel
#   90000  AGTCC>A   deletion, novel
#
# Usage: make_input.sh <out.vcf> [reference.fa]
set -euo pipefail
OUT="$1"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF="${2:-$HERE/../functional/fixtures/ref/tiny.fa}"

base_at() { samtools faidx "$REF" "chr20:$1-$2" | tail -n +2 | tr -d '\n' | tr '[:lower:]' '[:upper:]'; }

{
	# resources/vcf_header.txt is only a STUB -- it starts at ##FORMAT and carries no ##fileformat
	# line. RUFUS.interpret writes ##fileformat/##fileDate itself before appending the stub
	# (RUFUS.interpret.cpp:5303). Without that first line bcftools rejects the file outright with
	# "unknown file type", the stage reads zero variants, and the test exits 0 having proved nothing.
	echo '##fileformat=VCFv4.1'
	grep -v '^##contig' "$HERE/../../resources/vcf_header.txt"
	echo '##contig=<ID=chr20,length=200000>'
	printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSUBJECT\n'
	emit() { printf 'chr20\t%s\t.\t%s\t%s\t100\tPASS\tHD=30_30_30\tGT:KDP:KRO:KAO\t0/1:30:15:15\n' "$1" "$2" "$3"; }
	emit 20000 "$(base_at 20000 20002)" CTG
	emit 60000 "$(base_at 60000 60000)" A
	emit 70000 "$(base_at 70000 70002)" CAC
	emit 80000 "$(base_at 80000 80000)" "$(base_at 80000 80000)TTTT"
	emit 90000 "$(base_at 90000 90004)" "$(base_at 90000 90000)"
} > "$OUT"
echo "wrote $OUT ($(grep -vc '^#' "$OUT") records)" >&2
