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
# Baseline was re-frozen when #98 landed: the canonical output stopped being atomized and gained
# CO_ATOMS, and one real variant was recovered (chr20:2172408 TA>GG, whose second atom the control
# carries -- the old any-atom behaviour emitted T>G alone, asserting a haplotype the subject does not
# have). Any failure here now is a real regression.
set -uo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../lib.sh"
R="${RUFUS_REG_TEST_FILES:-$REPO_ROOT/../resources/reg_test_files}"
RUN="$R/runs/chr20_m5"
GEN=SMHTCOLO829T-X-X-M45-A001-uwsc-SMAFIW2W6WLI-sentieon_bwamem_202308.01_GRCh38.aligned.sorted.cram.chr20.generator
SUBJ=SMHTCOLO829T-X-X-M45-A001-uwsc-SMAFIW2W6WLI-sentieon_bwamem_202308.01_GRCh38.aligned.sorted.cram
CTRL="$R/crams/SMHTCOLO829BL-X-X-M45-A001-uwsc-SMAFIR3ZFBTN-sentieon_bwamem_202308.01_GRCh38.aligned.sorted.cram"
# The golden lives IN THE REPO, not in the run directory. Two reasons: re-freezing a baseline then
# shows up as a reviewable diff, which is the whole point of separating it from the no-op refactor;
# and the preserved run output stays untouched as the record of what the pre-#98 pipeline produced.
GOLD="$TESTS_DIR/golden/chr20_m5.canonical.vcf.gz"
GOLD_ATOMIZED="$TESTS_DIR/golden/chr20_m5.atomized.vcf.gz"
for f in "$RUN/rufus_chr20/$GEN.V2.overlap.hashcount.fastq.bam.vcf" "$CTRL" "$GOLD" "$R/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"; do
	[ -e "$f" ] || { echo "  SKIP: missing $f"; exit 0; }
done

TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
root="$TMP/r"; mkdir -p "$root/rufus_chr20/Intermediates"
# The preserved fixture predates the DP/RO/AO -> KDP/KRO/KAO rename, and that rename is a clean
# version boundary: the pipeline deliberately refuses a legacy VCF rather than guessing. So migrate
# the staged COPY to the current tag vocabulary -- the test then exercises the supported path, and
# reg_test_files stays pristine (it is the oracle; see #100 for what happens when it is not).
# This transformation is value-preserving: only tag names change, never numbers.
sed -e 's/GT:DP:RO:AO/GT:KDP:KRO:KAO/g' \
    -e 's/;AO=/;KAO=/g' \
    -e '/^##FORMAT=<ID=AK,/d' \
    -e 's/^##FORMAT=<ID=DP,\(.*\)$/##FORMAT=<ID=KDP,\1/' \
    -e 's/^##FORMAT=<ID=RO,\(.*\)$/##FORMAT=<ID=KRO,\1/' \
    -e 's/^##FORMAT=<ID=AO,\(.*\)$/##FORMAT=<ID=KAO,\1/' \
    -e 's/^##INFO=<ID=AO,\(.*\)$/##INFO=<ID=KAO,\1/' \
    "$RUN/rufus_chr20/$GEN.V2.overlap.hashcount.fastq.bam.vcf" \
    > "$root/rufus_chr20/$GEN.V2.overlap.hashcount.fastq.bam.vcf"

GEN="$GEN" SUBJ="$SUBJ" REF="$R/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa" \
EXTRA_BIND="$(cd "$R" && pwd)" \
	run_finalize "$root" --controls "$CTRL" >"$root/log" 2>&1
code=$?
[ "$code" = 77 ] && { echo "  SKIP: no runner"; exit 0; }
NEW="$root/temp.RUFUS.Final.$SUBJ.chr20.vcf.gz"
[ -s "$NEW" ] || { fail "no final VCF (exit $code); see $root/log"; exit 1; }

# Volatile header lines, stripped exactly as tests/replay/replay.sh:165 does:
#   ##bcftools_*        embed the command line and temp paths
#   ##fileDate          a timestamp from whenever RUFUS.interpret last produced the input
#   ##RUFUSCommandLine  likewise
# ##fileDate matters more than it looks: the input for this test is a preserved interpret output, and
# anything that regenerates it moves that line. The golden must not be sensitive to it.
volatile() { grep -v '^##bcftools' | grep -v '^##fileDate=' | grep -v '^##RUFUSCommandLine='; }
if ! diff <(bcftools view -H "$GOLD" 2>/dev/null) <(bcftools view -H "$NEW" 2>/dev/null) > "$TMP/d"; then
	fail "records differ from golden ($(wc -l < "$TMP/d") diff lines)"; head -6 "$TMP/d" >&2; exit 1
fi
if ! diff <(bcftools view -h "$GOLD" 2>/dev/null | volatile) \
          <(bcftools view -h "$NEW"  2>/dev/null | volatile) > "$TMP/dh"; then
	fail "headers differ beyond ##bcftools_* provenance"; head -6 "$TMP/dh" >&2; exit 1
fi
ok "canonical matches golden ($(bcftools view -H "$GOLD" 2>/dev/null | wc -l) records)"

# The atomized sidecar is part of the contract too (#98) -- consumers whose comparison is
# position-based point at it, so a silent change there is as bad as one in the canonical file.
NEW_ATOMIZED="$root/temp.RUFUS.Final.$SUBJ.chr20.atomized.vcf.gz"
[ -s "$NEW_ATOMIZED" ] || { fail "atomized sidecar was not produced"; exit 1; }
if ! diff <(bcftools view -H "$GOLD_ATOMIZED" 2>/dev/null) <(bcftools view -H "$NEW_ATOMIZED" 2>/dev/null) > "$TMP/da"; then
	fail "atomized sidecar differs from golden ($(wc -l < "$TMP/da") diff lines)"; head -6 "$TMP/da" >&2; exit 1
fi
ok "atomized sidecar matches golden ($(bcftools view -H "$GOLD_ATOMIZED" 2>/dev/null | wc -l) records)"
