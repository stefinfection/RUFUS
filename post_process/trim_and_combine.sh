#!/bin/bash

# This script concatenates individual vcfs produced by piecemeal RUFUS runs,
# into a single combined vcf. Only variants from within the region targeted by the run
# are kept (at this time, RUFUS may make erroneous calls out of the given range). 
# The final vcf contains the header from the first chr1_1_{CHUNK_SIZE} final 
# individual vcf. Note, the chr and chr length arrays, as well as the chunk size
# must be equal to those run for the piecemeal run.
# Compresses and indexes the final file.

module load bcftools
module load htslib

# Files and dirs for runscript
SOURCE_DIR=/scratch/ucgd/lustre-labs/marth/scratch/u0746015/HapMap/rufus_runs/yu/kmer_tests/whole_genome/k20_one_perc_test/

RUFUS_VCF=NA12878-NA18517-1percent-lib1-umi-liquid_tumor.sorted.bam.generator.V2.overlap.hashcount.fastq.bam.coinherited.vcf.gz
TRIMMED_VCF=trimmed.final.vcf
INHERITED_VCF=NA12878-NA18517-1percent-lib1-umi-liquid_tumor.sorted.bam.generator.FINAL.normalized.vcf.gz

RUFUS_PRE_VCF=Intermediates/NA12878-NA18517-1percent-lib1-umi-liquid_tumor.sorted.bam.generator.V2.overlap.hashcount.fastq.bam.sorted.vcf
TRIMMED_PRE_VCF=Intermediates/trimmed.prefiltered.vcf

RUFUS_SAMPLE_STRING="NA12878-NA18517-1percent-lib1-umi-liquid_tumor\tNA12878"

COMBINED_VCF="test.combined.final.vcf"
COMBINED_PRE_VCF="test.combined.prefiltered.vcf"

# Headers that get written to vcf
COMBINED_HEADER="combined.header"
COMBINED_PRE_HEADER="combined.preheader"

# Start of headers
HEADER_START="/scratch/ucgd/lustre-labs/marth/scratch/u0746015/HapMap/resources/combined.header.start"
PRE_HEADER_START="/scratch/ucgd/lustre-labs/marth/scratch/u0746015/HapMap/resources/combined.preheader.start"

# Records that get written to vcf (non-header)
COMBINED_RECORDS="combined.records"
COMBINED_PRE_RECORDS="combined.prerecords"

NUM_CHRS=24
CHUNK_SIZE=1000000

CHRS=(
"1"  
"2" 
"3" 
"4" 
"5" 
"6" 
"7" 
"8" 
"9" 
"10" 
"11" 
"12" 
"13" 
"14" 
"15" 
"16" 
"17" 
"18" 
"19" 
"20" 
"21" 
"22" 
"X" 
"Y" 
)

# GRCh38
CHR_LENGTHS=( 
248956422  
242193529 
198295559 
190214555 
181538259 
170805979 
159345973 
145138636 
138394717 
133797422 
135086622 
133275309 
114364328 
107043718 
101991189 
90338345 
83257441 
80373285 
58617616 
64444167 
46709983 
50818468 
156040895 
57227415 
)

cd $SOURCE_DIR

# Initialize combined headers
cat "$HEADER_START" > $COMBINED_HEADER
cat "$PRE_HEADER_START" > $COMBINED_PRE_HEADER
TEMP_TRIMMED="temp.trimmed"

n=${#CHRS[@]}
for (( i = 0; i < n; i++ ))
do
    curr_len="${CHR_LENGTHS[i]}"
    curr_chr="${CHRS[i]}"
    ((chunks=curr_len/CHUNK_SIZE+1))
 
    start_coord=1
    for (( c = 1; c <= chunks; c++ ))
    do
        ((end_coord=start_coord+CHUNK_SIZE-1))
        if [ $end_coord -gt $curr_len ]; then
            echo "passed end coordinate, setting to $curr_len"
            end_coord=$curr_len
        fi
	
	    # Get region dir
        TARGET_DIR="chr${curr_chr}/chr${curr_chr}_${start_coord}_${end_coord}/"

        if [[ -f "$TARGET_DIR${RUFUS_VCF}" ]]; then
       
            # Write out trimmed region to final vcf
            bcftools view -r "chr${curr_chr}:${start_coord}-${end_coord}" "${TARGET_DIR}${RUFUS_VCF}" > $TEMP_TRIMMED
            bcftools view -H $TEMP_TRIMMED >> $COMBINED_RECORDS
            bcftools view -h $TEMP_TRIMMED | grep "##contig" >> $COMBINED_HEADER 
            
            # Write out trimmed region to prefiltered vcf
            if [[ -f "${TARGET_DIR}${RUFUS_PRE_VCF}" ]]; then
                bgzip "${TARGET_DIR}${RUFUS_PRE_VCF}"
                bcftools index -t "${TARGET_DIR}${RUFUS_PRE_VCF}.gz"
            fi
            bcftools view -r "chr${curr_chr}:${start_coord}-${end_coord}" "${TARGET_DIR}${RUFUS_PRE_VCF}.gz" > $TEMP_TRIMMED
            bcftools view -H $TEMP_TRIMMED >> $COMBINED_PRE_RECORDS
            bcftools view -h $TEMP_TRIMMED | grep "##contig" >> $COMBINED_PRE_HEADER 
        fi
    
        # Advance start coordinate
        ((start_coord=end_coord+1))
    done
    echo "Done with chr${curr_chr}"
done

cat $COMBINED_HEADER | uniq > $COMBINED_VCF
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$RUFUS_SAMPLE_STRING" >> $COMBINED_VCF
cat $COMBINED_RECORDS >> $COMBINED_VCF

cat $COMBINED_PRE_HEADER | uniq > $COMBINED_PRE_VCF
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$RUFUS_SAMPLE_STRING" >> $COMBINED_PRE_VCF
cat $COMBINED_PRE_RECORDS >> $COMBINED_PRE_VCF

bgzip $COMBINED_VCF
bcftools index -t "${COMBINED_VCF}.gz"

bgzip $COMBINED_PRE_VCF
bcftools index -t "${COMBINED_PRE_VCF}.gz"

# Clean up temp files
rm $TEMP_TRIMMED
rm $COMBINED_HEADER
rm $COMBINED_PRE_HEADER
rm $COMBINED_RECORDS
rm $COMBINED_PRE_RECORDS
