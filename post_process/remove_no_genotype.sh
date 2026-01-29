#!/bin/bash
# Works for any number of sample columns in VCF
input_file=$1

# Process the gzipped VCF file
zcat "$input_file" | awk -F'\t' '
BEGIN {
    skipped_lines = 0;
    expected_cols = 0;
}
{
    # If this is a header line
    if ($0 ~ /^#/) {
        print
        # If this is the column header line, count the expected number of columns
        if ($0 ~ /^#CHROM/) {
            expected_cols = NF;
        }
    }
    # If this is a data line
    else {
        # Check if we have the expected number of columns (all genotype fields present)
        if (NF == expected_cols) {
            print
        } else {
            skipped_lines++
        }
    }
}
END {
    # Optional: report skipped lines to stderr
    if (skipped_lines > 0) {
        print "Skipped lines:", skipped_lines > "/dev/stderr"
    }
}'