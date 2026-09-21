#!/bin/bash
# A provider is conformant when it reproduces the golden table AND satisfies the contract rules.
# Both halves matter: the golden catches drift, the rules catch a provider that is self-consistently
# wrong.  Usage: test_conformance.sh [provider-name]
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
FIX="$HERE/fixture"; F="$ROOT/tests/functional/fixtures"
PROVIDER="${1:-bcftools}"
GOLD="$FIX/golden.subject.tsv"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
rc=0

bash "$ROOT/post_process/pileup/run_pileup.sh" --provider "$PROVIDER" \
	--sites-vcf "$FIX/sites.vcf" --bam "$F/somatic/tumor.bam" --ref "$F/ref/tiny.fa" \
	--role SUBJECT --sample-name tumor --out "$TMP/out.tsv" 2>"$TMP/err" || {
		echo "  FAIL: run_pileup.sh exited non-zero"; sed 's/^/    /' "$TMP/err"; exit 1; }

# The #provider and #params preamble lines legitimately differ between engines and settings.
strip() { grep -v '^#provider=\|^#params=' "$1"; }
if diff <(strip "$GOLD") <(strip "$TMP/out.tsv") > "$TMP/d"; then
	echo "  ok: reproduces the golden table ($(grep -vc '^#\|^CHROM' "$GOLD") sites)"
else
	echo "  FAIL: differs from golden ($(grep -c '^[<>]' "$TMP/d") lines)"; head -6 "$TMP/d" | sed 's/^/    /'; rc=1
fi

# Contract rule 4: the annotator must never meet a value it cannot parse. Every field is a number
# or "."; a provider inventing "NA" or "-" would break consumers in ways the golden diff would show
# only for this one fixture.
bad=$(grep -v '^#\|^CHROM' "$TMP/out.tsv" | awk -F'\t' '{for(i=6;i<=NF;i++) if($i !~ /^(\.|-?[0-9]+(\.[0-9]+)?(e[-+]?[0-9]+)?)$/) {print; break}}' | head -3)
if [ -n "$bad" ]; then echo "  FAIL: non-numeric, non-'.' values present"; echo "$bad" | sed 's/^/    /'; rc=1
else echo "  ok: every field is numeric or '.'"; fi

# The documented blind spots, asserted so an improvement is NOTICED rather than absorbed silently.
# See CONTRACT.md "Known limits". If these start passing, the contract needs updating -- that is a
# good outcome, but it must not happen quietly.
for spec in "150000:composite MNV" "170000:1000bp deletion"; do
	pos="${spec%%:*}"; what="${spec#*:}"
	alt=$(awk -F'\t' -v p="$pos" '$2==p {print $8}' "$TMP/out.tsv")
	if [ "${alt:-0}" = "0" ]; then echo "  ok: $what at $pos still unscored (AD_ALT=0), as documented"
	else echo "  NOTE: $what at $pos now scores AD_ALT=$alt -- provider improved; update CONTRACT.md"; rc=1; fi
done
exit $rc
