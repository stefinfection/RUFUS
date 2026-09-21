#!/bin/bash
# Run the finalize_vcf.sh suite.  Usage: run_all.sh [case-name-substring]
#
# A case marked "# XFAIL:" in its header is expected to fail until the named change lands. When one
# starts passing the runner reports XPASS -- that is the signal to delete the marker, not to ignore it.
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
echo "pass=$pass fail=$fail xfail=$xfail xpass=$xpass"
[ "$fail" -eq 0 ]
