#!/bin/bash
# bcftools pileup provider. See ../CONTRACT.md.
#
# Emits canonical rows only; run_pileup.sh writes the preamble and column header.
#
# Every non-obvious flag here is load-bearing, and each one was established by testing against a real
# 107GB CRAM rather than from the documentation:
#
#   -C alleles -T FILE  constrains counting to the caller's REF/ALT instead of whatever alleles the
#                       engine observes. The targets file wants CHROM<TAB>POS<TAB>REF,ALT -- REF and
#                       ALT comma-joined in ONE column. Separate columns fail outright.
#   -i                  keeps sites present in -T that mpileup did not observe. Without it a site with
#                       no read support VANISHES instead of reporting zeros, which contract rule 1
#                       forbids: a missing row and a zero row mean different things.
#   bcftools norm       REQUIRED as the final stage. mpileup rewrites indels into full repeat context
#                       (A>AT comes back as ATTTTTT>ATTTTTTT); 18 of 49 indels differed from input
#                       until this restored the caller's representation. Contract rule 2.
#   -B                  disables BAQ, which otherwise zeroes base qualities near indels and suppresses
#                       support for exactly the variants RUFUS most cares about.
#   -a ...              only the NON-default tags are requested. The starred tags in `mpileup -a '?'`
#                       (BQBZ, MQBZ, MQSBZ, RPBZ, SCBZ, MQ0F, SGB, VDB) arrive automatically and
#                       naming them explicitly is an error: "Could not parse tag INFO/MQSBZ".
set -euo pipefail

BAM=""; REF=""; ALLELES=""; REGIONS=""; ROLE=""; DEPTH=10000; MIN_MQ=0; MIN_BQ=13; BAQ=off
while [ $# -gt 0 ]; do
	case "$1" in
		--version) echo "bcftools:$(bcftools --version | head -1 | awk '{print $2}')"; exit 0;;
		--bam) BAM="$2"; shift 2;;      --ref) REF="$2"; shift 2;;
		--alleles) ALLELES="$2"; shift 2;; --regions) REGIONS="$2"; shift 2;;
		--role) ROLE="$2"; shift 2;;    --depth) DEPTH="$2"; shift 2;;
		--min-mq) MIN_MQ="$2"; shift 2;; --min-bq) MIN_BQ="$2"; shift 2;;
		--baq) BAQ="$2"; shift 2;;
		*) echo "bcftools provider: unknown argument '$1'" >&2; exit 2;;
	esac
done
for v in BAM REF ALLELES REGIONS ROLE; do
	[ -n "${!v}" ] || { echo "bcftools provider: missing --${v,,}" >&2; exit 2; }
done

baq_flag=(-B); [ "$BAQ" = on ] && baq_flag=()

bcftools mpileup -f "$REF" -R "$REGIONS" \
	-a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,FORMAT/SCR,FORMAT/NMBZ \
	-d "$DEPTH" -q "$MIN_MQ" -Q "$MIN_BQ" "${baq_flag[@]}" "$BAM" -Ou 2>/dev/null \
  | bcftools call -m -A -C alleles -T "$ALLELES" -i -Ou 2>/dev/null \
  | bcftools norm -f "$REF" -m- -Ou 2>/dev/null \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%DP]\t[%AD]\t[%ADF]\t[%ADR]\t%MQ0F\t%RPBZ\t%BQBZ\t%MQBZ\t%MQSBZ\t%SCBZ\t%SGB\t[%SP]\t[%SCR]\t[%NMBZ]\n' 2>/dev/null \
  | awk -F'\t' -v OFS='\t' -v role="$ROLE" '
	# Split the Number=R arrays into their REF and ALT halves. A site bcftools emitted with no
	# data at all gives "." for the whole array; report that as 0 support, not as missing --
	# the reads were looked at and said nothing.
	function half(arr, i,   a, n) { n = split(arr, a, ","); return (n >= i && a[i] != "" && a[i] != ".") ? a[i] : 0 }
	function num(v) { return (v == "" || v == ".") ? "." : v }
	{
		print $1, $2, $3, $4, role,
		      ($5 == "." ? 0 : $5),
		      half($6,1), half($6,2),
		      half($7,1), half($7,2), half($8,1), half($8,2),
		      ".", ".", ".", ".",              # F1R2/F2R1: bcftools cannot produce these (contract rule 4)
		      num($9), num($10), num($11), num($12), num($13), num($14), num($15),
		      num($16), ($17 == "." ? 0 : $17), num($18)
	}'
