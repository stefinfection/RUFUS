#!/bin/bash
# The all-atoms rule (#98): drop a composite record only if EVERY atom matched a control variant.
# Partial inheritance of a haplotype does not make the haplotype inherited.
#
# Reproduces the fabricated-variant failure exactly, on fixture data, in under a second:
#
#   input                 chr20:20000  AAC>CTG          one MNP, one haplotype
#   control germline      chr20:20000  A>C              carried by fixtures/somatic/normal.bam
#   decompose + isec -w1  drops atom 1, keeps atoms 2,3
#   OUTPUT TODAY          chr20:20001 A>T + 20002 C>G   <- asserts reference A at 20000
#
# The subject's haplotype is CTG. The emitted records claim a haplotype where 20000 is reference,
# which the subject does not have. That is not a dropped variant, it is a manufactured one.
#
# Expected after #98: the composite survives whole, annotated CO_ATOMS=1/3.
set -uo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../lib.sh"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
GEN=t.generator
CONTROL="$FIXTURES/somatic/normal.bam"
[ -e "$CONTROL" ] || { echo "  SKIP: missing $CONTROL"; exit 0; }

root="$TMP/r"; mkdir -p "$root/rufus_chr20/Intermediates"
bash "$TESTS_DIR/make_input.sh" "$root/rufus_chr20/$GEN.V2.overlap.hashcount.fastq.bam.vcf" 2>/dev/null
run_finalize "$root" --controls "$CONTROL" >"$root/log" 2>&1
code=$?
[ "$code" = 77 ] && { echo "  SKIP: no runner"; exit 0; }
out="$root/temp.RUFUS.Final.t.bam.chr20.vcf.gz"
[ -s "$out" ] || { fail "no final VCF (exit $code)"; exit 1; }

rc=0
span=$(repr "$out" | awk -F'\t' '$2>=20000 && $2<=20002')
if echo "$span" | awk -F'\t' '$2==20000 && length($3)==3 && length($4)==3' | grep -q .; then
	ok "composite chr20:20000 AAC>CTG survives whole"
else
	fail "composite was split; emitted: $(echo "$span" | tr '\n' ' ')"; rc=1
fi
if echo "$span" | awk -F'\t' '$2==20001 || $2==20002' | grep -q .; then
	fail "orphan atoms emitted at 20001/20002 -- these assert a haplotype the subject does not carry"
	rc=1
else
	ok "no orphan atoms left behind"
fi
co=$(bcftools query -f '%POS\t%INFO/CO_ATOMS\n' "$out" 2>/dev/null | awk -F'\t' '$1==20000{print $2}')
if [ -n "$co" ] && [ "$co" != "." ]; then ok "CO_ATOMS=$co"; else echo "  note: CO_ATOMS not present yet"; fi
exit $rc
