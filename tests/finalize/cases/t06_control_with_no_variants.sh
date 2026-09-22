#!/bin/bash
# A control that calls NO variants at the subject's sites must leave every record intact.
#
# This is the regression for the failure that took out functional tests f1, f2, f3 and f7: the
# co-inheritance keep-list was built with awk's NR==FNR two-file idiom, which does the wrong thing
# when the FIRST file is empty -- and empty is the NORMAL case here, because a control with no
# variants at these sites matches no atoms. FNR restarts at 1 for file 2 while NR is also 1, so the
# second file's rows were swallowed as if they were the first: every parent looked fully matched,
# every call was dropped, the emitted VCF came out zero bytes, and bcftools refused to index it --
# failing the run with exit 255 rather than emitting anything.
#
# The existing cases could not catch it. t02 passes --controls "" so remove_coinheriteds never runs,
# and t04/t05 use a real control that DOES call variants. Nothing exercised "control ran and found
# nothing", which is exactly what small fixtures produce.
set -uo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../lib.sh"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
GEN=t.generator
CONTROL="$FIXTURES/somatic/normal.bam"
[ -e "$CONTROL" ] || { echo "  SKIP: missing $CONTROL"; exit 0; }

root="$TMP/r"; mkdir -p "$root/rufus_chr20/Intermediates"
IN="$root/rufus_chr20/$GEN.V2.overlap.hashcount.fastq.bam.vcf"
bash "$TESTS_DIR/make_input.sh" "$IN" 2>/dev/null
# Drop chr20:20000 -- it sits on the control's germline SNV, so keeping it would give the control
# something to call and defeat the whole point of this case.
grep -v '^chr20	20000	' "$IN" > "$IN.tmp" && mv "$IN.tmp" "$IN"
want=$(grep -vc '^#' "$IN")

run_finalize "$root" --controls "$CONTROL" >"$root/log" 2>&1
code=$?
[ "$code" = 77 ] && { echo "  SKIP: no runner"; exit 0; }

rc=0
if ! grep -q 'zero variant records' "$root/log"; then
	echo "  NOTE: the control called variants after all; this case is not exercising the empty path"
fi
out="$root/temp.RUFUS.Final.t.bam.chr20.vcf.gz"
if [ ! -s "$out" ]; then
	echo "  FAIL: no final VCF (exit $code). An empty keep-list must still yield a valid VCF."
	grep -i 'keeping\|cannot be usefully indexed' "$root/log" | head -3 | sed 's/^/    /'
	exit 1
fi
got=$(bcftools view -H "$out" 2>/dev/null | wc -l)
if [ "$got" -eq "$want" ]; then
	echo "  ok: control called nothing, all $got record(s) survived"
else
	echo "  FAIL: $want record(s) in, $got out -- a control that matched nothing dropped calls"
	grep -i 'Co-inheritance: keeping' "$root/log" | sed 's/^/    /'
	rc=1
fi
# The emitted VCF must be valid even when the keep-list IS legitimately empty.
if bcftools view -h "$out" >/dev/null 2>&1; then echo "  ok: output is a valid, indexable VCF"
else echo "  FAIL: output is not a readable VCF"; rc=1; fi
exit $rc
