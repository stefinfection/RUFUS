#!/bin/bash
# Works for any number of sample columns in VCF
input_file=$1
output_file=$2

# Process the gzipped VCF file
zcat "$input_file" | awk -F'\t' '
BEGIN {
    skipped_lines = 0;
    expected_cols = 0;
}
{
    # If this is a header line
    if ($0 ~ /^#/) {
        print $0 > output_file;
        # If this is the column header line, count the expected number of columns
        if ($0 ~ /^#CHROM/) {
            expected_cols = NF;
        }
    } 
    # If this is a data line
    else {
        # Check if we have the expected number of columns (all genotype fields present)
        if (NF == expected_cols) {
            print $0 > output_file;
        } else {
            skipped_lines++;
        }
    }
}
END {
    print skipped_lines " lines were not printed because they did not have all genotype columns.";
}
' output_file="$output_file"

bgzip "$output_file"