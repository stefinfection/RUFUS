#!/bin/bash
# Test-only provider. Wraps bcftools and then deliberately (a) drops several fields to "." and
# (b) perturbs one count, so compare_providers.sh can be shown to DETECT both -- a comparison tool
# that only ever agrees has not been tested.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
[ "${1:-}" = --version ] && { echo "degraded:test"; exit 0; }
"$HERE/../../../post_process/pileup/providers/bcftools.sh" "$@" \
  | awk -F'\t' -v OFS='\t' '{ $17="."; $18="."; $19="."; if ($2==120000) $8=$8+1; print }'
