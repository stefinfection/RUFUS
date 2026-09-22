#!/bin/bash
# Put canonical pileup tables into a VCF. See CONTRACT.md.
#
# This script knows nothing about how the counts were produced -- no pileup flags, no provider names,
# no engine assumptions. It reads the table and the role mapping in the table's own preamble. That is
# the point of the split: replacing the engine must not touch this file.
#
# PHASE 1 IS ANNOTATE ONLY. Nothing here filters, genotypes, or alters an existing value. Every tag is
# additive, so records are unchanged apart from new FORMAT fields.
set -euo pipefail

VCF=""; OUT=""; TABLES=()
while [ $# -gt 0 ]; do
	case "$1" in
		--vcf)   VCF="$2"; shift 2;;
		--out)   OUT="$2"; shift 2;;
		--table) TABLES+=("$2"); shift 2;;
		-h|--help) sed -n '2,9p' "$0" >&2; exit 2;;
		*) echo "annotate_from_pileup.sh: unknown argument '$1'" >&2; exit 2;;
	esac
done
[ -n "$VCF" ] && [ -n "$OUT" ] && [ ${#TABLES[@]} -gt 0 ] || {
	echo "usage: annotate_from_pileup.sh --vcf IN.vcf.gz --out OUT.vcf.gz --table T.tsv [--table ...]" >&2; exit 2; }

TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT

# Canonical table column numbers (CONTRACT.md). Number=R tags occupy a REF column and the ALT column
# after it, and are recombined here into the comma-joined form VCF wants.
declare -A COL=( [DP]=6 [AD]=7 [ADF]=9 [ADR]=11 [F1R2]=13 [F2R1]=15
                 [MQ0F]=17 [RPBZ]=18 [BQBZ]=19 [MQBZ]=20 [MQSBZ]=21 [SCBZ]=22 [SGB]=23
                 [SP]=24 [SCR]=25 [NMBZ]=26 )
PAIRED=" AD ADF ADR F1R2 F2R1 "
ORDER=(DP AD ADF ADR F1R2 F2R1 MQ0F RPBZ BQBZ MQBZ MQSBZ SCBZ SGB SP SCR NMBZ)

# DP/AD/ADF/ADR/AF carry their SPEC meanings here -- read counts, which is what the KDP/KRO/KAO rename
# freed these names for, and what lets hap.py / vcfeval / GATK consume this output. The bias
# statistics keep bcftools' names but are FORMAT rather than INFO, because the pileup runs one
# invocation per sample (CONTRACT.md) and these values are therefore per-sample.
declare -A HDR=(
 [DP]='##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth from the pileup stage (reads, not k-mers; the k-mer statistic is KDP)">'
 [AD]='##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele from the pileup stage">'
 [ADF]='##FORMAT=<ID=ADF,Number=R,Type=Integer,Description="Read depth for each allele on the forward strand">'
 [ADR]='##FORMAT=<ID=ADR,Number=R,Type=Integer,Description="Read depth for each allele on the reverse strand">'
 [AF]='##FORMAT=<ID=AF,Number=1,Type=Float,Description="Read-based allele fraction AD[alt]/(AD[ref]+AD[alt]). ADVISORY ONLY for composite alleles and deletions longer than a read: the pileup engine cannot score those and reports 0 support even when reads carry the variant (see CONTRACT.md)">'
 [F1R2]='##FORMAT=<ID=F1R2,Number=R,Type=Integer,Description="Reads in F1R2 orientation supporting each allele">'
 [F2R1]='##FORMAT=<ID=F2R1,Number=R,Type=Integer,Description="Reads in F2R1 orientation supporting each allele">'
 [MQ0F]='##FORMAT=<ID=MQ0F,Number=1,Type=Float,Description="Fraction of reads with mapping quality zero">'
 [RPBZ]='##FORMAT=<ID=RPBZ,Number=1,Type=Float,Description="Mann-Whitney U-z of read position bias">'
 [BQBZ]='##FORMAT=<ID=BQBZ,Number=1,Type=Float,Description="Mann-Whitney U-z of base quality bias">'
 [MQBZ]='##FORMAT=<ID=MQBZ,Number=1,Type=Float,Description="Mann-Whitney U-z of mapping quality bias">'
 [MQSBZ]='##FORMAT=<ID=MQSBZ,Number=1,Type=Float,Description="Mann-Whitney U-z of mapping quality vs strand bias">'
 [SCBZ]='##FORMAT=<ID=SCBZ,Number=1,Type=Float,Description="Mann-Whitney U-z of soft-clip length bias">'
 [SGB]='##FORMAT=<ID=SGB,Number=1,Type=Float,Description="Segregation-based metric">'
 [SP]='##FORMAT=<ID=SP,Number=1,Type=Integer,Description="Phred-scaled strand bias p-value">'
 [SCR]='##FORMAT=<ID=SCR,Number=1,Type=Integer,Description="Number of soft-clipped reads at this site">'
 [NMBZ]='##FORMAT=<ID=NMBZ,Number=1,Type=Float,Description="Mann-Whitney U-z of mismatch count in supporting reads">'
 [AD_OTHER]='##FORMAT=<ID=AD_OTHER,Number=1,Type=Integer,Description="Reads attributable to NEITHER listed allele: DP - AD[ref] - AD[alt]. Deliberately named for the measurement, not a cause -- a composite allele the engine cannot score is the known case (its variant-carrying reads land here), but an indel where a substitution was asked about, or an allele the record does not name, would look the same. AD_OTHER large alongside AD[alt]=0 is the fingerprint of the engine being blind rather than the variant being absent">'
)
# Allele-level, so INFO rather than FORMAT. Records are biallelic here (norm -m- ran upstream).
NO_PILEUP_MODEL_HDR='##INFO=<ID=NO_PILEUP_MODEL,Number=1,Type=Integer,Description="1 when the pileup engine has no model for this allele class, so its read counts are not evidence about this allele. Derived from the SHAPE of REF/ALT, not from any pileup output: currently set for composite (equal-length multi-base) alleles, which are counted per-position and therefore never attributed to the composite. AF is reported missing rather than 0 for these records. NOTE: deletions longer than a read are equally unmodelled but are NOT flagged, as detecting them needs a read-length threshold -- see CONTRACT.md">'


# Does column $2 of table $1 hold a real value at any site?
#
# Deliberately ONE awk reading the file directly, with no pipeline and no early exit. The obvious
# formulation -- pipe the rows into `awk '$c != "." {exit}'` -- is wrong under `set -o pipefail` and
# wrong ONLY AT SCALE: awk quits on the first match, the upstream writer gets SIGPIPE, pipefail turns
# that into a failed test, and every tag is judged absent. A 9-row fixture passes because upstream
# finishes before awk exits; a 1108-site run reports DP=138 as "never present" and silently drops
# every provider tag.
has_value() {
	awk -F'\t' -v c="$2" '!/^#/ && $1 != "CHROM" && $c != "." { found = 1 } END { exit !found }' "$1"
}

CUR="$VCF"
for T in "${TABLES[@]}"; do
	[ -s "$T" ] || { echo "annotate_from_pileup.sh: missing or empty table $T" >&2; exit 1; }
	ver=$(awk -F= '/^#contract_version=/{print $2; exit}' "$T")
	[ "$ver" = 1 ] || { echo "annotate_from_pileup.sh: $T is contract version '${ver:-none}', expected 1" >&2; exit 1; }

	# Role mapping comes from the table, never from parsing a sample name out of a filename -- the
	# naming-by-filename coupling is a known source of WG-vs-region divergence elsewhere in RUFUS.
	ROLE=$(awk -F'\t' '/^#(SUBJECT|CONTROL_)/{sub(/^#/,"",$1); print $1; exit}' "$T")
	SNAME=$(awk -F'\t' '/^#(SUBJECT|CONTROL_)/{print $4; exit}' "$T")
	[ -n "$ROLE" ] || { echo "annotate_from_pileup.sh: $T has no '#<ROLE>' preamble line" >&2; exit 1; }
	if [ -z "$SNAME" ] || [ "$SNAME" = "." ]; then
		SNAME=$(bcftools query -l "$CUR" | head -1)   # SUBJECT is the first sample by RUFUS convention
	fi
	bcftools query -l "$CUR" | grep -qxF "$SNAME" || {
		echo "annotate_from_pileup.sh: sample '$SNAME' (role $ROLE) is not present in $CUR" >&2; exit 1; }

	# Contract rule 4: a provider may emit "." for anything it cannot compute. A tag that is "." at
	# EVERY site is dropped entirely rather than declared and left empty -- a header line for a field
	# nothing carries is exactly the dead declaration the KDP rename removed elsewhere. The bcftools
	# provider cannot produce F1R2/F2R1, so this path runs on every real invocation rather than
	# sitting untested until another engine needs it.
	present=(); exprs=(); spec="CHROM,POS,REF,ALT"
	for tag in "${ORDER[@]}"; do
		c=${COL[$tag]}
		if has_value "$T" "$c"; then
			present+=("$tag"); spec="$spec,FORMAT/$tag"
			case "$PAIRED" in *" $tag "*) exprs+=("\$$c\",\"\$$((c+1))");; *) exprs+=("\$$c");; esac
		fi
	done
	# AF is derived here rather than by the provider: it is a ratio of two numbers the provider already
	# reported, and computing it in both places invites the two to disagree.
	present+=(AF); spec="$spec,FORMAT/AF"; exprs+=("af")
	# Two additions that cost nothing and remove a trap: AF=0 for an unscorable allele is
	# indistinguishable from AF=0 for a variant that genuinely is not there.
	present+=(AD_OTHER); spec="$spec,FORMAT/AD_OTHER"
	spec="$spec,INFO/NO_PILEUP_MODEL"

	{ echo -n 'BEGIN{FS=OFS="\t"} /^#/||$1=="CHROM"{next} {'
	  echo -n 'r=($7=="."?0:$7); a=($8=="."?0:$8); t=r+a;'
	  echo -n 'd=($6=="."?0:$6); unacc=d-t; if(unacc<0) unacc=0;'
	  # Composite allele: REF and ALT the same length and longer than one base. Structural, so no
	  # threshold and no guessing -- either it is a multi-base substitution or it is not.
	  echo -n 'nomodel=(length($3)==length($4) && length($3)>1) ? 1 : 0;'
	  # AF is OUR derived number. Computing 0 from a numerator we know the engine could not fill is
	  # manufacturing a misleading value -- "." is the VCF-native way to say no information, and every
	  # consumer already skips missing without having to know about the flag.
	  echo -n 'af = nomodel ? "." : ((t>0) ? sprintf("%.4f",a/t) : ".");'
	  echo -n 'print $1,$2,$3,$4'
	  for e in "${exprs[@]}"; do echo -n ",$e"; done
	  echo ',unacc,nomodel}'
	} > "$TMP/build.awk"

	awk -f "$TMP/build.awk" "$T" | sort -k1,1 -k2,2n | bgzip > "$TMP/annot.$ROLE.tsv.gz"
	tabix -s1 -b2 -e2 -f "$TMP/annot.$ROLE.tsv.gz"
	: > "$TMP/hdr.$ROLE.txt"
	for tag in "${present[@]}"; do printf '%s\n' "${HDR[$tag]}" >> "$TMP/hdr.$ROLE.txt"; done
	printf '%s\n' "$NO_PILEUP_MODEL_HDR" >> "$TMP/hdr.$ROLE.txt"

	NEXT="$TMP/annotated.$ROLE.vcf.gz"
	bcftools annotate -s "$SNAME" -a "$TMP/annot.$ROLE.tsv.gz" -h "$TMP/hdr.$ROLE.txt" \
		-c "$spec" "$CUR" -Oz -o "$NEXT"
	bcftools index -f "$NEXT"
	echo "annotate_from_pileup: $ROLE -> sample '$SNAME', ${#present[@]} tags: ${present[*]}" >&2
	CUR="$NEXT"
done

cp "$CUR" "$OUT"
bcftools index -f "$OUT"
