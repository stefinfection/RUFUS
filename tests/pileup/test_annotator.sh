#!/bin/bash
# The annotator must add tags and change nothing else. Phase 1 is annotate-only, so "the records are
# otherwise untouched" is the property that makes it safe to land before any filtering exists.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
FIX="$HERE/fixture"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
rc=0

# A sites VCF with a sample column to annotate into.
bcftools view "$FIX/sites.vcf" -Ov 2>/dev/null | awk -F'\t' -v OFS='\t' '
	/^##fileformat/{print; print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"; next}
	/^##/{print; next} /^#CHROM/{print $0,"FORMAT","tumor"; next} {print $0,"GT","0/1"}' \
	| bgzip > "$TMP/in.vcf.gz"
bcftools index -f "$TMP/in.vcf.gz"

bash "$ROOT/post_process/pileup/annotate_from_pileup.sh" --vcf "$TMP/in.vcf.gz" \
	--out "$TMP/out.vcf.gz" --table "$FIX/golden.subject.tsv" 2>"$TMP/err" \
	|| { echo "  FAIL: annotator exited non-zero"; sed 's/^/    /' "$TMP/err"; exit 1; }

# 1. ANNOTATE ONLY: every pre-existing column must be byte-identical. Compare columns 1-8, which is
#    everything except FORMAT and the sample -- the only places new tags may appear.
if diff <(bcftools view -H "$TMP/in.vcf.gz" 2>/dev/null | cut -f1-8) \
        <(bcftools view -H "$TMP/out.vcf.gz" 2>/dev/null | cut -f1-8 | sed 's/;PU_UNSCORED=[01]//; s/^\([^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t[^\t]*\t\)PU_UNSCORED=[01]$/\1./') >/dev/null; then
	echo "  ok: CHROM..INFO unchanged apart from the added PU_UNSCORED"
else
	echo "  FAIL: the annotator altered a pre-existing column"; rc=1
fi

# 2. Values reach the VCF intact.
# The table has a multi-line '#' preamble AND a column-header row; both must be skipped or the
# comparison picks up the literal string "POS  DP  AD_REF,AD_ALT" as if it were a site.
awk -F'\t' '!/^#/ && $1 != "CHROM" { print $2 "\t" $6 "\t" $7 "," $8 }' "$FIX/golden.subject.tsv" \
	| sort > "$TMP/want"
bcftools query -f '%POS\t[%DP\t%AD]\n' "$TMP/out.vcf.gz" 2>/dev/null | sort > "$TMP/got"
if diff "$TMP/want" "$TMP/got" > "$TMP/vd"; then
	echo "  ok: DP and AD match the canonical table at every site ($(wc -l < "$TMP/want") sites)"
else
	echo "  FAIL: values disagree with the table"; head -4 "$TMP/vd" | sed 's/^/    /'; rc=1
fi

# 3. Contract rule 4: a tag the provider never computed must not be declared. bcftools cannot produce
#    F1R2/F2R1, so this exercises the degradation path on a real run rather than a contrived one.
n=$(bcftools view -h "$TMP/out.vcf.gz" 2>/dev/null | grep -c 'ID=F1R2\|ID=F2R1')
if [ "$n" -eq 0 ]; then echo "  ok: all-'.' tags dropped rather than declared empty"
else echo "  FAIL: $n header line(s) declare a tag the provider never produced"; rc=1; fi

# 4. The unscorable-allele marker. AF=0 on a composite allele means "the engine cannot answer", and
#    without this flag it is indistinguishable from "no reads support it".
u=$(bcftools query -f '%POS\t%INFO/PU_UNSCORED\n' "$TMP/out.vcf.gz" 2>/dev/null | awk -F'\t' '$2==1 {print $1}' | tr '\n' ' ')
if [ "$u" = "150000 " ]; then echo "  ok: PU_UNSCORED set on the composite allele only"
else echo "  FAIL: PU_UNSCORED expected on 150000 only, got: ${u:-none}"; rc=1; fi

# 5. AF is a ratio of the reported counts, not an independent estimate.
af=$(bcftools query -f '%POS\t[%AF]\n' "$TMP/out.vcf.gz" 2>/dev/null | awk -F'\t' '$1==120000{print $2}')
if [ "$af" = "0.4722" ]; then echo "  ok: AF derived from AD (17/36 = 0.4722)"
else echo "  FAIL: AF at 120000 expected 0.4722, got '$af'"; rc=1; fi
exit $rc
