#!/bin/bash
# Qualify one provider against another over the same input.
#
# This is how a replacement engine gets trusted with real data: not by reading its code, but by
# showing it agrees with the incumbent site by site. Counts must match exactly -- a read either
# supports an allele or it does not. The bias statistics are floating point and may legitimately
# differ in the last digits, so those get a tolerance.
#
# Usage: compare_providers.sh <providerA> <providerB> [tolerance]
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
FIX="$HERE/fixture"; F="$ROOT/tests/functional/fixtures"
A="${1:?provider A}"; B="${2:?provider B}"; TOL="${3:-0.001}"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT

for p in "$A" "$B"; do
	bash "$ROOT/post_process/pileup/run_pileup.sh" --provider "$p" \
		--sites-vcf "$FIX/sites.vcf" --bam "$F/somatic/tumor.bam" --ref "$F/ref/tiny.fa" \
		--role SUBJECT --out "$TMP/$p.tsv" 2>"$TMP/$p.err" \
		|| { echo "FAIL: provider '$p' did not run"; sed 's/^/  /' "$TMP/$p.err"; exit 1; }
done

awk -F'\t' -v tol="$TOL" -v na="$A" -v nb="$B" '
	function isnum(v) { return v ~ /^-?[0-9]+(\.[0-9]+)?(e[-+]?[0-9]+)?$/ }
	NR==FNR { if (/^#/) next; if ($1=="CHROM") { for(i=1;i<=NF;i++) h[i]=$i; next } a[$1":"$2":"$3":"$4]=$0; next }
	/^#/ { next }
	$1=="CHROM" { next }
	{
		k=$1":"$2":"$3":"$4
		if (!(k in a)) { printf "  site only in %s: %s\n", nb, k; bad++; next }
		split(a[k], av, "\t"); seen[k]=1
		for (i=6; i<=NF; i++) {
			x=av[i]; y=$i
			if (x==y) continue
			# "." vs a value is a CAPABILITY difference, not a disagreement: one engine computes
			# the field and the other does not. Report it separately -- it is expected during a swap.
			if (x=="." || y==".") { cap[h[i]]++; continue }
			if (isnum(x) && isnum(y)) {
				d = x-y; if (d<0) d=-d
				# counts (columns 6-16) must match exactly; float stats get the tolerance
				if (i<=16) { printf "  %s  %s: %s=%s %s=%s (counts must match exactly)\n", k, h[i], na, x, nb, y; bad++ }
				else if (d>tol) { printf "  %s  %s: %s=%s %s=%s (delta %.4g > %s)\n", k, h[i], na, x, nb, y, d, tol; bad++ }
			} else { printf "  %s  %s: %s=%s %s=%s\n", k, h[i], na, x, nb, y; bad++ }
		}
	}
	END {
		for (k in a) if (!(k in seen)) { printf "  site only in %s: %s\n", na, k; bad++ }
		if (length(cap)) { printf "\n  capability differences (one side emits \".\"):\n"; for (c in cap) printf "    %-8s %d site(s)\n", c, cap[c] }
		printf "\n  %s\n", bad ? "DISAGREE: " bad " difference(s)" : "AGREE on every site"
		exit bad ? 1 : 0
	}' "$TMP/$A.tsv" "$TMP/$B.tsv"
