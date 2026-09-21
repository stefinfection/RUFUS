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
 [PU_UNACC]='##FORMAT=<ID=PU_UNACC,Number=1,Type=Integer,Description="Reads at this site supporting NEITHER listed allele (DP - AD[ref] - AD[alt]). A large value means the reads carry something the engine could not attribute to REF or ALT -- for a composite allele that is the variant itself">'
)
# Allele-level, so INFO rather than FORMAT. Records are biallelic here (norm -m- ran upstream).
PU_UNSCORED_HDR='##INFO=<ID=PU_UNSCORED,Number=1,Type=Integer,Description="1 when the pileup engine cannot score this allele class, so AD/AF are meaningless for it rather than merely zero. Currently set for composite (equal-length multi-base) alleles, which mpileup has no model for. NOTE: deletions longer than a read are also unscorable but are NOT flagged here, because that needs a read-length threshold -- see CONTRACT.md">'


rows() { grep -v '^#' "$1" | awk -F'\t' 'NR>1 || $1!="CHROM"' | awk -F'\t' '$1!="CHROM"'; }

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
		if rows "$T" | awk -F'\t' -v c="$c" '$c != "." {found=1; exit} END{exit !found}'; then
			present+=("$tag"); spec="$spec,FORMAT/$tag"
			case "$PAIRED" in *" $tag "*) exprs+=("\$$c\",\"\$$((c+1))");; *) exprs+=("\$$c");; esac
		fi
	done
	# AF is derived here rather than by the provider: it is a ratio of two numbers the provider already
	# reported, and computing it in both places invites the two to disagree.
	present+=(AF); spec="$spec,FORMAT/AF"; exprs+=("af")
	# Two additions that cost nothing and remove a trap: AF=0 for an unscorable allele is
	# indistinguishable from AF=0 for a variant that genuinely is not there.
	present+=(PU_UNACC); spec="$spec,FORMAT/PU_UNACC"
	spec="$spec,INFO/PU_UNSCORED"

	{ echo -n 'BEGIN{FS=OFS="\t"} /^#/||$1=="CHROM"{next} {'
	  echo -n 'r=($7=="."?0:$7); a=($8=="."?0:$8); t=r+a; af=(t>0)?sprintf("%.4f",a/t):".";'
	  echo -n 'd=($6=="."?0:$6); unacc=d-t; if(unacc<0) unacc=0;'
	  # composite allele: REF and ALT the same length and longer than one base. Structural, so no
	  # threshold and no guessing -- either it is a multi-base substitution or it is not.
	  echo -n 'unscored=(length($3)==length($4) && length($3)>1) ? 1 : 0;'
	  echo -n 'print $1,$2,$3,$4'
	  for e in "${exprs[@]}"; do echo -n ",$e"; done
	  echo ',unacc,unscored}'
	} > "$TMP/build.awk"

	awk -f "$TMP/build.awk" "$T" | sort -k1,1 -k2,2n | bgzip > "$TMP/annot.$ROLE.tsv.gz"
	tabix -s1 -b2 -e2 -f "$TMP/annot.$ROLE.tsv.gz"
	: > "$TMP/hdr.$ROLE.txt"
	for tag in "${present[@]}"; do printf '%s\n' "${HDR[$tag]}" >> "$TMP/hdr.$ROLE.txt"; done
	printf '%s\n' "$PU_UNSCORED_HDR" >> "$TMP/hdr.$ROLE.txt"

	NEXT="$TMP/annotated.$ROLE.vcf.gz"
	bcftools annotate -s "$SNAME" -a "$TMP/annot.$ROLE.tsv.gz" -h "$TMP/hdr.$ROLE.txt" \
		-c "$spec" "$CUR" -Oz -o "$NEXT"
	bcftools index -f "$NEXT"
	echo "annotate_from_pileup: $ROLE -> sample '$SNAME', ${#present[@]} tags: ${present[*]}" >&2
	CUR="$NEXT"
done

cp "$CUR" "$OUT"
bcftools index -f "$OUT"
