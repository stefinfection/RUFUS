#!/bin/bash
# Golden replay: the stage must reproduce a preserved real run byte for byte.
#
# This is the regression test for the extraction itself and for every future change to the stage -- it
# is the thing that makes "no-op refactor" a falsifiable claim rather than an assertion. It exercises
# the remove_coinheriteds branch against a real control CRAM, which the tiny fixtures cannot.
#
# LOCAL TIER: needs resources/reg_test_files (a sibling of the repo, not version-controlled) and the
# ~113GB control CRAM. Skips cleanly when they are absent.
#
# NOTE: after #98 lands, the canonical output stops being atomized, so this baseline MUST be
# re-frozen as part of that change. A failure here after #98 is expected exactly once; a failure at
# any other time is a real regression.
set -uo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../lib.sh"
R="${RUFUS_REG_TEST_FILES:-$REPO_ROOT/../resources/reg_test_files}"
RUN="$R/runs/chr20_m5"
GEN=SMHTCOLO829T-X-X-M45-A001-uwsc-SMAFIW2W6WLI-sentieon_bwamem_202308.01_GRCh38.aligned.sorted.cram.chr20.generator
SUBJ=SMHTCOLO829T-X-X-M45-A001-uwsc-SMAFIW2W6WLI-sentieon_bwamem_202308.01_GRCh38.aligned.sorted.cram
CTRL="$R/crams/SMHTCOLO829BL-X-X-M45-A001-uwsc-SMAFIR3ZFBTN-sentieon_bwamem_202308.01_GRCh38.aligned.sorted.cram"
GOLD="$RUN/temp.RUFUS.Final.$SUBJ.chr20.vcf.gz"
for f in "$RUN/rufus_chr20/$GEN.V2.overlap.hashcount.fastq.bam.vcf" "$CTRL" "$GOLD" "$R/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"; do
	[ -e "$f" ] || { echo "  SKIP: missing $f"; exit 0; }
done

TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
root="$TMP/r"; mkdir -p "$root/rufus_chr20/Intermediates"
cp "$RUN/rufus_chr20/$GEN.V2.overlap.hashcount.fastq.bam.vcf" "$root/rufus_chr20/"

GEN="$GEN" SUBJ="$SUBJ" REF="$R/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa" \
EXTRA_BIND="$(cd "$R" && pwd)" \
	run_finalize "$root" --controls "$CTRL" >"$root/log" 2>&1
code=$?
[ "$code" = 77 ] && { echo "  SKIP: no runner"; exit 0; }
NEW="$root/temp.RUFUS.Final.$SUBJ.chr20.vcf.gz"
[ -s "$NEW" ] || { fail "no final VCF (exit $code); see $root/log"; exit 1; }

# ##bcftools_* header lines embed the command line and temp paths, so they differ every run by design.
if ! diff <(bcftools view -H "$GOLD" 2>/dev/null) <(bcftools view -H "$NEW" 2>/dev/null) > "$TMP/d"; then
	fail "records differ from golden ($(wc -l < "$TMP/d") diff lines)"; head -6 "$TMP/d" >&2; exit 1
fi
if ! diff <(bcftools view -h "$GOLD" 2>/dev/null | grep -v '^##bcftools') \
          <(bcftools view -h "$NEW"  2>/dev/null | grep -v '^##bcftools') > "$TMP/dh"; then
	fail "headers differ beyond ##bcftools_* provenance"; head -6 "$TMP/dh" >&2; exit 1
fi
ok "reproduces chr20_m5 golden exactly ($(bcftools view -H "$GOLD" 2>/dev/null | wc -l) records)"
