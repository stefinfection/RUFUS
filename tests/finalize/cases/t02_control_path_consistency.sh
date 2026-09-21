#!/bin/bash
# THE issue #98 guard: the representation RUFUS emits must not depend on how controls were supplied.
#
# A BAM/CRAM control runs remove_coinheriteds.sh, which applies `vt normalize | vt decompose_blocksub`
# and whose output reaches the final VCF through `isec -w1`. A hash/generator/fastq control skips that
# entirely (runRufus.sh:1810). Today `bcftools norm -a` also decomposes MNVs, so the two paths
# converge and this test passes. Remove atomize without handling decompose_blocksub and they diverge:
# block substitutions stay split on the BAM path and intact on the other. That is the regression this
# test exists to catch.
#
# The two runs legitimately call DIFFERENT variant sets -- only the BAM path removes co-inherited
# variants -- so the assertion is on the representation of the variants present in BOTH, never on set
# equality.
set -uo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../lib.sh"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
GEN=t.generator
CONTROL="$FIXTURES/somatic/normal.bam"
[ -e "$CONTROL" ] || { echo "  SKIP: missing $CONTROL"; exit 0; }

for arm in withbam nocontrol; do
	root="$TMP/$arm"; mkdir -p "$root/rufus_chr20/Intermediates"
	bash "$TESTS_DIR/make_input.sh" "$root/rufus_chr20/$GEN.V2.overlap.hashcount.fastq.bam.vcf" 2>/dev/null
	if [ "$arm" = withbam ]; then ctrl="$CONTROL"; else ctrl=""; fi
	run_finalize "$root" --controls "$ctrl" >"$root/log" 2>&1
	code=$?
	[ "$code" = 77 ] && { echo "  SKIP: no runner"; exit 0; }
	out="$root/temp.RUFUS.Final.t.bam.chr20.vcf.gz"
	[ -s "$out" ] || { fail "$arm produced no final VCF (exit $code); see $root/log"; cp -r "$root" /tmp/t02_"$arm"_debug 2>/dev/null; exit 1; }
	repr "$out" > "$TMP/$arm.repr"
done

# Compare representation only where both arms kept the same site.
join -t"$(printf '\t')" -j1 \
	<(awk -F'\t' '{print $1":"$2"\t"$3"\t"$4}' "$TMP/withbam.repr" | sort) \
	<(awk -F'\t' '{print $1":"$2"\t"$3"\t"$4}' "$TMP/nocontrol.repr" | sort) > "$TMP/shared"
shared=$(wc -l < "$TMP/shared")
[ "$shared" -gt 0 ] || { fail "no sites shared between the two arms -- test proves nothing"; exit 1; }
bad=$(awk -F'\t' '$2!=$4 || $3!=$5' "$TMP/shared")
if [ -n "$bad" ]; then
	fail "representation differs by control type at $(echo "$bad" | wc -l) of $shared shared site(s):"
	echo "$bad" | head -5 | awk -F'\t' '{print "        "$1"  bam-control: "$2">"$3"   hash-control: "$4">"$5}' >&2
	exit 1
fi
ok "representation identical across control types at all $shared shared site(s)"
