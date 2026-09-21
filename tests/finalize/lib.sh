# Shared helpers for the finalize_vcf.sh test suite.
#
# These tests exercise post_process/finalize_vcf.sh directly rather than through a full RUFUS run.
# The stage takes ~30s standalone against real data and under a second against the tiny fixtures,
# versus hours for an end-to-end run, which is the whole reason it was extracted (issue #98).

TESTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$TESTS_DIR/../.." && pwd)"
FIXTURES="$REPO_ROOT/tests/functional/fixtures"
FINALIZE="$REPO_ROOT/post_process/finalize_vcf.sh"

# remove_coinheriteds.sh needs `vt`, which is only in the container. Run natively when vt is on PATH,
# otherwise wrap in singularity. Override the image with RUFUS_SIF.
RUFUS_SIF="${RUFUS_SIF:-/scratch/ucgd/lustre-labs/marth/resources/software/RUFUS/zenodo_images/rufus_latest.sif}"

# run_finalize <work_root> <extra finalize args...>
# Creates <work_root>/rufus_<region>/ and runs the stage against an input VCF the caller has already
# placed there as <generator>.V2.overlap.hashcount.fastq.bam.vcf.
run_finalize() {
	local work_root="$1"; shift
	local cmd=(bash "$FINALIZE" --work-dir "$work_root/rufus_${FR:-chr20}" --work-root "$work_root"
	           --generator "${GEN:-t.generator}" --subject-name "${SUBJ:-t.bam}"
	           --formatted-region "${FR:-chr20}" --region-postfix "${POSTFIX:-.chr20}"
	           --ref "${REF:-$FIXTURES/ref/tiny.fa}" --rufus-root "$REPO_ROOT"
	           --region "${REGION:-chr20}" --mosaic "${MOSAIC:-TRUE}" --threads "${THREADS:-2}"
	           "$@")
	if command -v vt >/dev/null 2>&1; then
		"${cmd[@]}"
	else
		[ -e "$RUFUS_SIF" ] || { echo "SKIP: no vt on PATH and no image at $RUFUS_SIF" >&2; return 77; }
		singularity exec --bind "$REPO_ROOT" --bind "$(dirname "$work_root")" \
			${EXTRA_BIND:+--bind $EXTRA_BIND} "$RUFUS_SIF" "${cmd[@]}"
	fi
}

# Representation of the variants in a VCF, as a sorted CHROM POS REF ALT table.
repr() { bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$1" 2>/dev/null | sort -k1,1 -k2,2n; }

fail() { echo "  FAIL: $*" >&2; return 1; }
ok()   { echo "  ok: $*"; }
