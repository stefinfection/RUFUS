#!/bin/bash
# finalize_vcf.sh must stop cleanly, with the right reason, on each "nothing to call" condition.
# The reason string is what runRufus.sh writes into region_status.log, so a wrong one here means a
# region is mislabelled in the run log -- silently, since these are all successful exits.
set -uo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../lib.sh"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
GEN=t.generator; rc=0

check() { # <name> <expected_reason> <expected_code> <setup>
	local name="$1" want_reason="$2" want_code="$3" setup="$4"
	local root="$TMP/$name"; mkdir -p "$root/rufus_chr20/Intermediates"
	local vcf="$root/rufus_chr20/$GEN.V2.overlap.hashcount.fastq.bam.vcf"
	eval "$setup"
	run_finalize "$root" --controls "" >"$root/log" 2>&1
	local got_code=$?
	[ "$got_code" = 77 ] && { echo "  SKIP: $name (no runner)"; return 0; }
	local got_reason; got_reason="$(cat "$root/rufus_chr20/.finalize_reason" 2>/dev/null || echo '<none>')"
	if [ "$got_reason" = "$want_reason" ] && [ "$got_code" = "$want_code" ]; then
		ok "$name -> $got_reason (exit $got_code)"
	else
		fail "$name: want '$want_reason'/$want_code, got '$got_reason'/$got_code"; rc=1
	fi
}

check missing_input no_interpret_passing_vars_dne 0 'true'
check empty_input   no_interpret_passing_vars_e   0 \
	'printf "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" > "$vcf"'
check unsanitizable no_records_after_sanitization 0 \
	'printf "##fileformat=VCFv4.2\n##contig=<ID=chr20,length=200000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr20\t100\t.\tA\t.\t.\t.\t.\n" > "$vcf"'
exit $rc
