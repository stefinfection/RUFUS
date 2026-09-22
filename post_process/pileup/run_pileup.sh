#!/bin/bash
# Produce a canonical pileup table for ONE sample. See CONTRACT.md.
#
# This script owns the contract: the preamble, the column header, and the site list. Providers own
# only "given reads and sites, what do the reads say" -- so swapping the engine cannot change the
# format, and a new provider cannot accidentally redefine a column.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONTRACT_VERSION=1

PROVIDER=bcftools; VCF=""; BAM=""; REF=""; ROLE=""; KIND=""; SAMPLE_NAME=""
DEPTH=10000; MIN_MQ=0; MIN_BQ=13; BAQ=off; OUT=""
while [ $# -gt 0 ]; do
	case "$1" in
		--provider) PROVIDER="$2"; shift 2;; --sites-vcf) VCF="$2"; shift 2;;
		--bam) BAM="$2"; shift 2;;           --ref) REF="$2"; shift 2;;
		--role) ROLE="$2"; shift 2;;         --kind) KIND="$2"; shift 2;;
		--sample-name) SAMPLE_NAME="$2"; shift 2;;
		--depth) DEPTH="$2"; shift 2;;       --min-mq) MIN_MQ="$2"; shift 2;;
		--min-bq) MIN_BQ="$2"; shift 2;;     --baq) BAQ="$2"; shift 2;;
		--out) OUT="$2"; shift 2;;
		-h|--help) sed -n '2,12p' "$0" >&2; exit 2;;
		*) echo "run_pileup.sh: unknown argument '$1'" >&2; exit 2;;
	esac
done
for v in VCF BAM REF ROLE; do
	[ -n "${!v}" ] || { echo "run_pileup.sh: missing --${v,,}" >&2; exit 2; }
done
# Tests supply their own providers (a deliberately degraded one, for instance) without adding them
# to the production provider directory.
# Search the test-supplied directory first, then the production one, so a test can add a provider
# (a deliberately degraded one, say) and still compare it against the real bcftools provider without
# copying or symlinking that into place.
PROVIDER_SH=""
for d in ${RUFUS_PILEUP_PROVIDER_DIR:+"$RUFUS_PILEUP_PROVIDER_DIR"} "$HERE/providers"; do
	[ -x "$d/$PROVIDER.sh" ] && { PROVIDER_SH="$d/$PROVIDER.sh"; break; }
done
[ -n "$PROVIDER_SH" ] || { echo "run_pileup.sh: no provider '$PROVIDER' in ${RUFUS_PILEUP_PROVIDER_DIR:+$RUFUS_PILEUP_PROVIDER_DIR or }$HERE/providers" >&2; exit 2; }
[ -n "$KIND" ] || KIND="$([ "${BAM##*.}" = cram ] && echo cram || echo bam)"

TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT

# Site list. `*` is the spanning-deletion placeholder: it denotes ABSENCE of sequence, so there is
# nothing to count and it is excluded here rather than handed to a provider that would have to
# invent an answer. The annotator reports those sites with a reason code (CONTRACT.md).
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$VCF" \
  | awk -F'\t' '$4 != "*"' | sort -k1,1 -k2,2n > "$TMP/sites"

# Two derived files, and they MUST come from the same source: with `call -i` the targets file drives
# which rows appear, so a site in one and not the other yields a phantom row with no data.
#   alleles: CHROM POS REF,ALT  -- the comma-joined third column `-C alleles` requires
#   regions: CHROM POS          -- what mpileup actually reads
awk -F'\t' -v OFS='\t' '{print $1, $2, $3","$4}' "$TMP/sites" | bgzip > "$TMP/alleles.tsv.gz"
tabix -s1 -b2 -e2 -f "$TMP/alleles.tsv.gz"
awk -F'\t' -v OFS='\t' '{print $1, $2}' "$TMP/sites" | sort -k1,1 -k2,2n -u > "$TMP/regions"

emit() {
	printf '#contract_version=%s\n' "$CONTRACT_VERSION"
	printf '#provider=%s\n' "$("$PROVIDER_SH" --version)"
	printf '#params=depth=%s,minMQ=%s,minBQ=%s,baq=%s\n' "$DEPTH" "$MIN_MQ" "$MIN_BQ" "$BAQ"
	printf '#role\tfile\tkind\tvcf_sample_name\n'
	printf '#%s\t%s\t%s\t%s\n' "$ROLE" "$BAM" "$KIND" "${SAMPLE_NAME:-.}"
	printf 'CHROM\tPOS\tREF\tALT\tROLE\tDP\tAD_REF\tAD_ALT\tADF_REF\tADF_ALT\tADR_REF\tADR_ALT'
	printf '\tF1R2_REF\tF1R2_ALT\tF2R1_REF\tF2R1_ALT\tMQ0F\tRPBZ\tBQBZ\tMQBZ\tMQSBZ\tSCBZ\tSGB\tSP\tSCR\tNMBZ\n'
	"$PROVIDER_SH" --bam "$BAM" --ref "$REF" --alleles "$TMP/alleles.tsv.gz" \
		--regions "$TMP/regions" --role "$ROLE" --depth "$DEPTH" \
		--min-mq "$MIN_MQ" --min-bq "$MIN_BQ" --baq "$BAQ"
}

if [ -n "$OUT" ]; then emit > "$OUT"; else emit; fi

# Contract rules 1 and 2, checked rather than trusted: both failure modes are silent. A provider that
# drops unsupported sites still produces a plausible table, and one that relabels a site still
# produces the right ROW COUNT -- so compare the actual keys, not just how many there are.
#
# Relabeling is not hypothetical. The trailing `bcftools norm` that restores mpileup's repeat-context
# rewriting also LEFT-ALIGNS, so a site supplied non-left-aligned comes back at a different POS
# (chr20:40000 ACAGTTTATTA>A returns as chr20:39999 CACAGTTTATT>C). Sites must already be normalized;
# this says so instead of silently producing a table the annotator cannot join.
if [ -n "$OUT" ]; then
	grep -v '^#' "$OUT" | tail -n +2 | cut -f1-4 | sort -k1,1 -k2,2n > "$TMP/out.keys"
	if ! diff -q "$TMP/sites" "$TMP/out.keys" >/dev/null 2>&1; then
		echo "run_pileup.sh: provider '$PROVIDER' did not return the sites it was given." >&2
		echo "  Supplied but missing from the output:" >&2
		comm -23 "$TMP/sites" "$TMP/out.keys" | head -5 | sed 's/^/    /' >&2
		echo "  Returned but never supplied (relabeled?):" >&2
		comm -13 "$TMP/sites" "$TMP/out.keys" | head -5 | sed 's/^/    /' >&2
		echo "  Sites must be left-aligned and trimmed before being passed here (CONTRACT.md rule 2)." >&2
		exit 1
	fi
fi
