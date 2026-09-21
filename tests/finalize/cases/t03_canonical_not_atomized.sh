#!/bin/bash
# XFAIL: until #98 stops atomizing the canonical output
# XFAIL until #98: the canonical VCF must keep composite alleles, not decompose them into atoms.
#
# Atomization is lossy -- it destroys the linkage a single contig asserted, inflates variant counts,
# and manufactures records no read supports (see #98 for the chr20:62459362 case, where an atom came
# back AD=136,0 at DP=141). RUFUS is assembly-based, so a composite allele is the caller being
# faithful to the haplotype it assembled.
#
# Planted: chr20:70000 ATG>CAC, a 3bp block substitution. Expected in the canonical output as ONE
# record. Today `bcftools norm -a` splits it into three SNVs, so this fails by design until #98 lands.
set -uo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../lib.sh"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
GEN=t.generator
root="$TMP/r"; mkdir -p "$root/rufus_chr20/Intermediates"
bash "$TESTS_DIR/make_input.sh" "$root/rufus_chr20/$GEN.V2.overlap.hashcount.fastq.bam.vcf" 2>/dev/null
run_finalize "$root" --controls "" >"$root/log" 2>&1
code=$?
[ "$code" = 77 ] && { echo "  SKIP: no runner"; exit 0; }
out="$root/temp.RUFUS.Final.t.bam.chr20.vcf.gz"
[ -s "$out" ] || { fail "no final VCF (exit $code)"; exit 1; }

rc=0
# 1. the block substitution survives whole
if repr "$out" | awk -F'\t' '$2==70000 && length($3)==3 && length($4)==3' | grep -q .; then
	ok "chr20:70000 kept as one composite record"
else
	fail "chr20:70000 ATG>CAC was decomposed; output has: $(repr "$out" | awk -F'\t' '$2>=70000 && $2<=70002' | tr '\n' ' ')"
	rc=1
fi
# 2. no star alleles -- with atomize gone, `*` should only ever come from genuinely overlapping
#    deletions between separate events, never as bookkeeping for our own decomposition
if repr "$out" | awk -F'\t' '$4 ~ /\*/' | grep -q .; then
	fail "output contains '*' alleles: $(repr "$out" | awk -F'\t' '$4 ~ /\*/' | head -3 | tr '\n' ' ')"
	rc=1
else
	ok "no '*' spanning-deletion placeholders in canonical output"
fi
exit $rc
