#!/bin/bash
# Finalize a RUFUS region/whole-genome VCF: everything between RUFUS.interpret's raw output and the
# per-invocation final VCF -- sanitize, dedupe, PASS gate, reference fix-up, region trim,
# co-inherited removal, normalization, HD_AF, sort.
#
# Extracted verbatim from runRufus.sh (the block that ran inline after the overlap stage) so that the
# representation-changing steps live in one callable place and can be tested without a full RUFUS run.
# See issue #98.
#
# SHELL SEMANTICS: runRufus.sh runs under `set -e` only (line 28) -- no `-u`, no `pipefail`. This
# script must match exactly, or extracting the block changes its behaviour. Do NOT "tidy" this to
# `set -euo pipefail`.
set -e

usage() {
	cat >&2 <<'EOF'
Usage: finalize_vcf.sh --work-dir DIR --work-root DIR --generator NAME --subject-name NAME \
                       --formatted-region STR --region-postfix STR --ref FASTA --rufus-root DIR \
                       [--region REGION] [--mosaic TRUE|FALSE] [--threads N] [--controls a,b,c]

Consumes  $WORK_DIR/$GENERATOR.V2.overlap.hashcount.fastq.bam.vcf
Produces  $WORK_ROOT/temp.RUFUS.Final.$SUBJECT_NAME$REGION_POSTFIX.vcf.gz (+ .csi)

Early-stop conditions write the reason to $WORK_DIR/.finalize_reason and exit with the code the
inline block used (0 for "nothing to call", 100 for tabix failures). The caller maps that back onto
its own _region_exit_reason so region_status.log is unchanged.
EOF
	exit 2
}

WORK_DIR=""; WORK_ROOT=""; ProbandGenerator=""; ProbandFileName=""
formatted_region=""; region_postfix=""; _arg_ref=""; RDIR=""
_arg_region=""; _arg_mosaic="FALSE"; _arg_threads=10; CONTROLS=()
PILEUP=auto; PILEUP_PROVIDER=bcftools; PILEUP_DEPTH=10000; SUBJECT_ALIGN=""

while [ $# -gt 0 ]; do
	case "$1" in
		--work-dir)         WORK_DIR="$2"; shift 2;;
		--work-root)        WORK_ROOT="$2"; shift 2;;
		--generator)        ProbandGenerator="$2"; shift 2;;
		--subject-name)     ProbandFileName="$2"; shift 2;;
		--formatted-region) formatted_region="$2"; shift 2;;
		--region-postfix)   region_postfix="$2"; shift 2;;
		--ref)              _arg_ref="$2"; shift 2;;
		--rufus-root)       RDIR="$2"; shift 2;;
		--region)           _arg_region="$2"; shift 2;;
		--mosaic)           _arg_mosaic="$2"; shift 2;;
		--threads)          _arg_threads="$2"; shift 2;;
		--controls)         if [ -n "$2" ]; then IFS=',' read -r -a CONTROLS <<< "$2"; fi; shift 2;;
		--pileup)           PILEUP="$2"; shift 2;;
		--pileup-provider)  PILEUP_PROVIDER="$2"; shift 2;;
		--pileup-depth)     PILEUP_DEPTH="$2"; shift 2;;
		--subject-align)    SUBJECT_ALIGN="$2"; shift 2;;
		-h|--help)          usage;;
		*) echo "finalize_vcf.sh: unknown argument '$1'" >&2; usage;;
	esac
done

for _req in WORK_DIR WORK_ROOT ProbandGenerator ProbandFileName formatted_region _arg_ref RDIR; do
	if [ -z "${!_req}" ]; then
		echo "finalize_vcf.sh: missing required argument for $_req" >&2
		usage
	fi
done

# Early-stop channel back to the caller. Cleared on entry so a stale file from a previous run in the
# same WORK_DIR can never be mistaken for this invocation's outcome.
# The sub-scripts this stage calls (remove_coinheriteds.sh, add_hd_med.add_hd_af.sh,
# remove_no_genotype.sh, single_pileup.sh) all read WORK_DIR from the ENVIRONMENT, not from an
# argument. runRufus.sh exports it; inheriting that only works when we are its child, so export it
# here as well or a standalone invocation dies inside the first sub-script.
export WORK_DIR

REASON_FILE="$WORK_DIR/.finalize_reason"
rm -f "$REASON_FILE"

stop_with() {
	printf '%s' "$1" > "$REASON_FILE"
	exit "$2"
}

# ----------------------------------------------------------------------------------------------
# Below this line the body is the runRufus.sh block verbatim, except that the five
# `_region_exit_reason=...; exit N` pairs became `stop_with ... N`.
# ----------------------------------------------------------------------------------------------

intermed_vcf="${WORK_DIR}/${ProbandGenerator}.V2.overlap.hashcount.fastq.bam.vcf"

if [[ -s "$intermed_vcf" ]]; then
	count=$(bcftools view -H "$intermed_vcf" | wc -l)
	if [ "$count" -eq 0 ]; then
	  	echo "Intermediate vcf contains no variants, indicating no variants found for this region." >&2
		stop_with "no_interpret_passing_vars_e" 0
	fi
  	# safe to proceed (file exists, non-empty, has variants)
else
	echo "Intermediate vcf not present, indicating no variants found for this region." >&2
	stop_with "no_interpret_passing_vars_dne" 0
fi

# Sanitize intermediate VCF: remove malformed records from RUFUS.Interpret output
# Keeps header lines as-is. For data lines, requires:
#   - at least 8 tab-delimited fields (CHROM POS ID REF ALT QUAL FILTER INFO)
#   - POS (col 2) is a positive integer
#   - REF (col 4) is non-empty and contains only valid bases (ACGTN)
#   - ALT (col 5) is non-empty, not just ".", and contains only valid VCF ALT characters
#   - If INFO (col 8) contains END=<n>, then END >= POS (prevents tabix "end < begin" error)
sanitized_vcf="${intermed_vcf}.sanitized.vcf"
awk -F'\t' '
/^#/ { print; next }
{
	if (NF < 8) next
	if ($2 !~ /^[0-9]+$/ || $2+0 < 1) next
	if ($4 == "" || $4 == "." || $4 !~ /^[ACGTNacgtn]+$/) next
	if ($5 == "" || $5 == ".") next
	if ($5 !~ /^[ACGTNacgtn.,*<>\[\]0-9:]+$/) next
	# Check END >= POS if END tag is present in INFO field
	info = $8
	if (match(info, /END=[0-9]+/)) {
		end_val = substr(info, RSTART+4, RLENGTH-4) + 0
		if (end_val < $2+0) next
	}
	print
}
' "$intermed_vcf" > "$sanitized_vcf"

sanitized_count=$(grep -vc "^#" "$sanitized_vcf" || true)
original_count=$(grep -vc "^#" "$intermed_vcf" || true)
removed_count=$((original_count - sanitized_count))
if [ "$removed_count" -gt 0 ]; then
	echo "WARNING: Removed $removed_count malformed VCF record(s) from RUFUS.Interpret output ($sanitized_count of $original_count records kept)." >&2
fi
if [ "$sanitized_count" -eq 0 ]; then
	echo "No valid VCF records remain after sanitization for this region." >&2
	stop_with "no_records_after_sanitization" 0
fi
mv "$sanitized_vcf" "$intermed_vcf"

# Trim off generator postfix
DEDUPED_VCF="$WORK_DIR/deduped.${formatted_region}.vcf"

# TODO: do I really need this? can I just sort?
grep "^#" "$intermed_vcf" > $WORK_DIR/Intermediates/$ProbandGenerator.V2.overlap.hashcount.fastq.bam.sorted.vcf
grep -v "^#" "$intermed_vcf" | sort -k1,1V -k2,2n >> $WORK_DIR/Intermediates/$ProbandGenerator.V2.overlap.hashcount.fastq.bam.sorted.vcf

echo "arg_mosaic = $_arg_mosaic"
if [ "$_arg_mosaic" == "TRUE" ]
then
	echo "including mosaic"
	bash $RDIR/scripts/VilterAutosomeOnly $WORK_DIR/Intermediates/$ProbandGenerator.V2.overlap.hashcount.fastq.bam.sorted.vcf "$_arg_ref" | perl $RDIR/scripts/ColapsDuplicateCalls.stream.pl > $DEDUPED_VCF
else
	echo "excluding mosaic"
	bash $RDIR/scripts/VilterAutosomeOnly.withoutMosaic $WORK_DIR/Intermediates/$ProbandGenerator.V2.overlap.hashcount.fastq.bam.sorted.vcf "$_arg_ref" | perl $RDIR/scripts/ColapsDuplicateCalls.stream.pl > $DEDUPED_VCF
fi

bgzip -f "$DEDUPED_VCF"
# Index with tabix, iteratively removing records that cause indexing failures
# This catches any malformed records that slip past the awk sanitizer
tabix_max_retries=50
tabix_attempt=0
while true; do
	tabix_stderr=$(tabix -C "$DEDUPED_VCF.gz" 2>&1) && break

	tabix_attempt=$((tabix_attempt + 1))
	if [ "$tabix_attempt" -ge "$tabix_max_retries" ]; then
		echo "ERROR: tabix failed after removing $tabix_attempt malformed record(s). Giving up." >&2
		echo "Last tabix error: $tabix_stderr" >&2
		stop_with "tabix_max_retries" 100
	fi

	# Parse the 1-based sequence number from: "Invalid record on sequence #N"
	bad_seq=$(echo "$tabix_stderr" | grep -oP 'sequence #\K[0-9]+' | head -1)
	if [ -z "$bad_seq" ]; then
		echo "ERROR: tabix failed with unexpected error: $tabix_stderr" >&2
		stop_with "tabix_unexpected_error" 100
	fi

	echo "WARNING: tabix indexing failed on data record #${bad_seq}, removing it and retrying (attempt $tabix_attempt)." >&2
	echo "  tabix error: $tabix_stderr" >&2

	# Decompress, remove the offending data line, recompress
	tmp_fix_vcf="${DEDUPED_VCF}.tabixfix.vcf"
	zcat "$DEDUPED_VCF.gz" | awk -v bad="$bad_seq" '
		/^#/ { print; next }
		{ data_line++; if (data_line != bad) print }
	' > "$tmp_fix_vcf"
	bgzip -f "$tmp_fix_vcf"
	mv "$tmp_fix_vcf.gz" "$DEDUPED_VCF.gz"
done

# Update reference alleles
REF_VCF="$WORK_DIR/ref.${formatted_region}.vcf"
bcftools +fill-from-fasta "$DEDUPED_VCF.gz" -- -c REF -f "$_arg_ref" > "$REF_VCF"

# Get rid of break-ends
TYPE_VCF="$WORK_DIR/snv_indel.${formatted_region}.vcf"
bcftools view -e "TYPE='bnd'" "$REF_VCF" > "$TYPE_VCF"

# Check for empty gt field
GX_VCF="$WORK_DIR/gx.${formatted_region}.vcf"
bash $RDIR/post_process/remove_no_genotype.sh "$TYPE_VCF" > "$GX_VCF"
bgzip -f "$GX_VCF"
bcftools index -f "$GX_VCF.gz"

# Trim calls to region. In whole-genome mode _arg_region is empty; `bcftools view -r ""` segfaults,
# and there is nothing to trim to, so pass the calls through unchanged.
TRIMMED_VCF="$WORK_DIR/trimed.${formatted_region}.vcf.gz"
if [ -n "$_arg_region" ]; then
	bcftools view -r "$_arg_region" "$GX_VCF.gz" -Oz -o "$TRIMMED_VCF"
else
	cp "$GX_VCF.gz" "$TRIMMED_VCF"
fi
bcftools index "$TRIMMED_VCF"

NO_CO_VCF="$WORK_DIR/no_coinheriteds.vcf.gz"
# remove_coinheriteds pileups each control at the variant sites, so it needs an alignable BAM/CRAM.
# A control given as a pre-built hash (a .generator stub) or fastq has no BAM to pile up (bwa would
# align an empty file -> mpileup fails on the empty bam). Collect only the BAM/CRAM controls and run
# the filter over those; skip entirely if none -- the HashList subtraction has already removed those
# controls' k-mers, so the co-inherited pileup is a secondary check with nothing to pile up.
_bamcram_controls=()
for _ctrl in "${CONTROLS[@]}"; do
	case "$_ctrl" in
		*.bam|*.cram) _bamcram_controls+=("$_ctrl") ;;
	esac
done
if [ ${#_bamcram_controls[@]} -ne "0" ]; then
	bash ${RDIR}/post_process/remove_coinheriteds.sh -t $_arg_threads -r "$formatted_region" -f "$_arg_ref" -i "$TRIMMED_VCF" -o "$NO_CO_VCF" -w "1000" -c "$(IFS=','; echo "${_bamcram_controls[*]}")"
else
	[ ${#CONTROLS[@]} -ne "0" ] && echo "Skipping remove_coinheriteds: no BAM/CRAM control to pile up (controls are hash/generator/fastq); HashList subtraction already handled them." >&2
	mv "$TRIMMED_VCF" "$NO_CO_VCF"
fi

# Left-align and split multiallelics. Deliberately NOT atomized (issue #98): atomization destroys the
# linkage a single contig asserted, inflates variant counts, and manufactures records no read
# supports. RUFUS is assembly-based, so a composite allele is the caller being faithful to the
# haplotype it actually assembled. An atomized copy is still emitted as a sidecar below, for
# consumers whose comparison is position/allele-string based rather than haplotype-aware.
CANON_VCF="$WORK_DIR/normalized.${formatted_region}.vcf"
bcftools norm -m- -f "$_arg_ref" "$NO_CO_VCF" -Oz -o "$CANON_VCF"

# Add HD_AF field
CANON_VCF_BASENAME=$(basename "$CANON_VCF")
HDAF_VCF="$WORK_DIR/hd_af.$CANON_VCF_BASENAME"
SUBJECT_SAMPLE_NAME=$(bcftools view -h "$CANON_VCF" | tail -n 1 | awk -F'\t' '{ print $10 }')
bash ${RDIR}/post_process/add_hd_med.add_hd_af.sh "$CANON_VCF" "$SUBJECT_SAMPLE_NAME" "$formatted_region"

# ---------------------------------------------------------------------------------------------
# Pileup annotation (phase 1 of the pileup work, #97). INFORMANT, NOT ADJUDICATOR: this adds read
# counts alongside the k-mer statistics and changes no FILTER, no genotype, and no existing value.
# Nothing downstream acts on these numbers yet -- that is deliberate, so the two channels can be
# compared on real data before either is trusted to decide anything.
#
# Subject only. Control pileups belong to phase 3, and leaving them out here also avoids having to
# map control roles onto VCF sample columns, which is reconstructed by filename surgery upstream and
# is a known source of WG-versus-region divergence.
# ---------------------------------------------------------------------------------------------
PILEUP_IN="$HDAF_VCF"
if [ "$PILEUP" != "off" ]; then
	# Resolve what to pile up. A bam/cram subject gives real reference reads; a generator or fastq
	# subject has none, so fall back to the aligned mutant reads, which carry ALT support only.
	_pu_align=""; _pu_kind=""
	if [ -n "$SUBJECT_ALIGN" ] && [ -e "$SUBJECT_ALIGN" ]; then
		_pu_align="$SUBJECT_ALIGN"
		case "$SUBJECT_ALIGN" in *.cram) _pu_kind=cram;; *) _pu_kind=bam;; esac
	elif [ -e "$WORK_DIR/$ProbandGenerator.Mutations.fastq.bam" ]; then
		_pu_align="$WORK_DIR/$ProbandGenerator.Mutations.fastq.bam"
		_pu_kind=altonly
		echo "Pileup: no aligned subject supplied; using the mutant-read BAM. ALT counts only -- it" >&2
		echo "        contains no reference reads, so AD[ref] and AF are not meaningful." >&2
	fi

	if [ -z "$_pu_align" ]; then
		echo "Pileup: skipped, no alignable subject input found." >&2
	else
		_pu_tsv="$WORK_DIR/pileup.SUBJECT.${formatted_region}.tsv"
		_pu_out="$WORK_DIR/pileup.${formatted_region}.vcf.gz"
		# A pileup failure must not lose the run: the k-mer VCF is complete without it, and phase 1
		# adds no value anything depends on. Warn and carry the un-annotated VCF forward.
		if bash "$RDIR/post_process/pileup/run_pileup.sh" \
				--provider "$PILEUP_PROVIDER" --sites-vcf "$PILEUP_IN" \
				--bam "$_pu_align" --ref "$_arg_ref" --role SUBJECT --kind "$_pu_kind" \
				--sample-name "$SUBJECT_SAMPLE_NAME" --depth "$PILEUP_DEPTH" \
				--out "$_pu_tsv" \
		   && bash "$RDIR/post_process/pileup/annotate_from_pileup.sh" \
				--vcf "$PILEUP_IN" --out "$_pu_out" --table "$_pu_tsv"; then
			PILEUP_IN="$_pu_out"
			echo "Pileup: annotated $(bcftools view -H "$_pu_out" 2>/dev/null | wc -l) record(s) from $_pu_kind input." >&2
		else
			echo "WARNING: pileup annotation failed; continuing without it. The k-mer calls are" >&2
			echo "         unaffected -- nothing downstream reads these tags yet." >&2
		fi
	fi
fi

# Sort
SORTED_VCF="$WORK_DIR/sorted.${formatted_region}.vcf.gz"
bcftools sort "$PILEUP_IN" -Oz -o "$SORTED_VCF"

# Rename final vcf and zip/index
PREFINAL_VCF="$WORK_DIR/temp.RUFUS.Final.${ProbandFileName}${region_postfix}.vcf.gz"
mv "$SORTED_VCF" "$PREFINAL_VCF"
bcftools index "$PREFINAL_VCF"

FINAL_BASENAME="$(basename "$PREFINAL_VCF")"

cp "$PREFINAL_VCF" "$WORK_ROOT/$FINAL_BASENAME"
cp "$PREFINAL_VCF.csi" "$WORK_ROOT/$FINAL_BASENAME.csi"

# Atomized sidecar (issue #98). The canonical VCF above keeps composite alleles; this copy decomposes
# them for downstream consumers whose comparison is position/allele-string based. Haplotype-aware
# comparison (GA4GH hap.py / RTG vcfeval) does not need it -- a composite allele and its decomposed
# equivalent compare EQUAL there -- so prefer the canonical file where the tooling allows.
# --old-rec-tag stamps each atom with a pointer back to the composite record it came from, so the
# sidecar is not a dead end.
ATOMIZED_BASENAME="temp.RUFUS.Final.${ProbandFileName}${region_postfix}.atomized.vcf.gz"
bcftools norm -a --old-rec-tag OLD_REC -f "$_arg_ref" "$PREFINAL_VCF" -Oz -o "$WORK_DIR/$ATOMIZED_BASENAME"
bcftools index -f "$WORK_DIR/$ATOMIZED_BASENAME"
cp "$WORK_DIR/$ATOMIZED_BASENAME" "$WORK_ROOT/$ATOMIZED_BASENAME"
cp "$WORK_DIR/$ATOMIZED_BASENAME.csi" "$WORK_ROOT/$ATOMIZED_BASENAME.csi"
