#!/bin/bash

# NOTE: THIS ONLY WORKS FOR A RUFUS RUN WITH TWO SAMPLE COLUMNS IN VCF - i.e. ONE NORMAL ONE TUMOR

input_file=$1
output_file=$2

# Process the gzipped VCF file
zcat "$input_file" | awk -F'\t' '
BEGIN {
    skipped_lines = 0;
}
{
    if ($0 ~ /^#/) {
        print $0 > output_file;
    } else if (NF == 11) {
        print $0 > output_file;
    } else {
        skipped_lines++;
    }
}
END {
    print skipped_lines " lines were not printed because they did not have the genotype columns.";
}
' output_file="$output_file"

bgzip $output_file
