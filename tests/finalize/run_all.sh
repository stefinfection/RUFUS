#!/bin/bash
# Run the finalize_vcf.sh suite.  Usage: run_all.sh [case-name-substring]
#
# Optional, and unused at present: a case marked "# XFAIL:" in its header is expected to fail until
# the named change lands, and the runner reports XPASS once it starts passing -- the signal to delete
# the marker. Used while building #98 to make the target behaviour executable before implementing it;
# kept for the next change that wants the same.
cd "$(dirname "${BASH_SOURCE[0]}")"
filter="${1:-}"
pass=0 fail=0 xfail=0 xpass=0
for c in cases/t*.sh; do
	name="$(basename "$c" .sh)"
	[ -n "$filter" ] && [[ "$name" != *"$filter"* ]] && continue
	xfail_reason="$(grep -m1 '^# XFAIL:' "$c" | sed 's/^# XFAIL: *//')"
	echo "== $name${xfail_reason:+  [XFAIL: $xfail_reason]}"
	if bash "$c"; then
		if [ -n "$xfail_reason" ]; then echo "  XPASS -- remove the XFAIL marker"; xpass=$((xpass+1))
		else pass=$((pass+1)); fi
	else
		if [ -n "$xfail_reason" ]; then echo "  (expected failure)"; xfail=$((xfail+1))
		else fail=$((fail+1)); fi
	fi
done
echo
summary="pass=$pass fail=$fail"
[ "$xfail" -gt 0 ] && summary="$summary xfail=$xfail"
[ "$xpass" -gt 0 ] && summary="$summary xpass=$xpass"
echo "$summary"
[ "$fail" -eq 0 ]
