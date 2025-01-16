#!/bin/bash

usage() {
	echo "Usage: $0 [-w window_size] [-r reference] [-s subject] [-c control1,control2,control3...] [-d source_dir] [-h]"
	echo "Options:"
	echo " -w window_size	Required: The size of the window used in the RUFUS run"
	echo " -r reference	Required: The reference used in the RUFUS run"
	echo " -c controls	Required: The control bam files used in the RUFUS run"
	echo " -s subject_file	Required: The name of the subject file: must be the same as that supplied to the RUFUS run"
	echo " -d source_dir	Required: The source directory where the RUFUS vcf(s) are located"
	echo " -h help	Print help message"
	exit 1
}

report_empty_results() {
  echo "RUFUS did not find any variants for the provided parameters. Please adjust and try again."
  echo "RUFUS did not find any variants for the provided parameters. Please adjust and try again." >&2
  echo "RUFUS did not find any variants for the provided parameters. Please adjust and try again." > results.out
}

# Cleans up intermediate files, reports no variants found in both out + error, and exits failure code
clean_up_early_intermeds() {
  local SUBJECT_FILE="$1"
  local ALL_ARGS=("$@")
  local CONTROLS=("${ALL_ARGS[@]:1}")

  # TODO: will need to not hard code eventually to accommodate other builds/species
  chroms=(
    "chr1",
    "chr2",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chr20",
    "chr21",
    "chr22",
    "chrX",
    "chrY",
  )

  # Clean up intermediate files
  echo -n "Cleaning up intermediates..." >&2

  # Have to do this piecemeal because too many files with windowed mode for single rm command
  for chrom in "${chroms[@]}"; do
    echo "looking for echo /mnt/${SUBJECT_FILE}*${chrom}*.generator*" >&2
    if ls /mnt/${SUBJECT_FILE}*${chrom}*.generator* 1> /dev/null 2>&1; then
      echo "found and trying to remove /mnt/${SUBJECT_FILE}*${chrom}*.generator*" >&2
      rm /mnt/${SUBJECT_FILE}*${chrom}*.generator*
    fi
  done

  for control in "${CONTROLS[@]}"; do
  echo -n "control $control" >&2
    # Have to do this piecemeal because too many files with windowed mode for single rm command
    for chrom in "${chroms[@]}"; do
        if ls /mnt/${control}*${chrom}*.generator* 1> /dev/null 2>&1; then
          echo "trying to remove /mnt/${control}*${chrom}*.generator*" >&2
          rm /mnt/${control}*${chrom}*.generator*
        fi
      done
  done

  # Remove intermediate files if they exist
  if [ -d "/mnt/Intermediates" ] && [ -d "/mnt/TempOverlap" ]; then
    rm -r /mnt/Intermediates
    rm -r /mnt/TempOverlap
  fi

  if [ -e "/mnt/temp*.vcf*" ]; then
    rm /mnt/temp*.vcf*
  fi
}

# static paths
bcftools="/opt/bcftools/bcftools"

# initialize vars
CONTROLS=()
WINDOW_SIZE=0
REFERENCE=""
SUBJECT_FILE=""
SOURCE_DIR="/mnt"

# parse command line arguments
while getopts ":w:r:c:s:d:h" option; do 
	case $option in 
		h) usage;;
		w) WINDOW_SIZE=$OPTARG;;
		r) REFERENCE=$OPTARG;;
		s) SUBJECT_FILE=$OPTARG;;
		d) SOURCE_DIR=$OPTARG;;
		c) IFS=',' read -r -a CONTROLS <<< "$OPTARG";;
		\?) echo "Invalid option: -$OPTARG" >&2
		    usage;;	
		*) echo "Option -$OPTARG requires an argument" >&2
		   usage;;
	esac
done
shift $((OPTIND-1))

# check for mandatory command line arguments
if [[ -z "$WINDOW_SIZE" ]]; then
	    echo "ERROR: Missing required option -w (window size)" >&2
fi

if [[ -z "$REFERENCE" ]]; then
	    echo "ERROR: Missing required option -r (reference)" >&2
fi

if [[ -z "$SOURCE_DIR" ]]; then
	echo "ERROR: Missing required option -d (source directory for RUFUS vcf(s))" >&2
fi

if [ ${#CONTROLS[@]} -eq 0 ]; then
	    echo "ERROR: Must supply at least one control bam" >&2
fi

cd $SOURCE_DIR
echo "RUFUS post-process version D-1.0.1"
date
start_time=$(date +"%s")

POST_PROCESS_DIR=/opt/RUFUS/post_process/
TEMP_FINAL_VCF="temp.RUFUS.Final.${SUBJECT_FILE}.combined.vcf.gz"
TEMP_PREFILTERED_VCF="temp.RUFUS.Prefiltered.${SUBJECT_FILE}.combined.vcf.gz"
GERMLINE_VCF="with_germline.RUFUS.Final.${SUBJECT_FILE}.combined.vcf.gz"

# Slight name change if not doing a windowed run
if [ "$WINDOW_SIZE" = "0" ]; then
	TEMP_FINAL_VCF="temp.RUFUS.Final.${SUBJECT_FILE}.vcf.gz"
	TEMP_PREFILTERED_VCF="${SUPP_DIR}temp.RUFUS.Prefiltered.${SUBJECT_FILE}.vcf.gz"
fi

# Check to see if final vcf exists, if not report empty results and exit
if [ ! -e "$TEMP_FINAL_VCF" ]; then
  report_empty_results
  clean_up_early_intermeds "$SUBJECT_FILE" "${CONTROLS[@]}"
fi

# Get number of variants reported
VARS_REPORTED=$($bcftools view -H $TEMP_FINAL_VCF | wc -l)

# Keep germline vcf
cp $TEMP_FINAL_VCF $GERMLINE_VCF
mv $GERMLINE_VCF rufus_supplementals/

# If windowed mode, trim and combine region
if [ "$WINDOW_SIZE" != "0" ]; then
	IFS=$'\t'
	TAB_DELIM_CONTROL_STRING="${CONTROLS[*]}"
	echo "Windowed run performed, trimming and combining region vcfs..."
	bash ${POST_PROCESS_DIR}trim_and_combine.sh $SUBJECT_FILE $TAB_DELIM_CONTROL_STRING $WINDOW_SIZE
fi

# Check for empty vcf AFTER trimming and combining
# If we don't have any variants here, the entire run didn't find any variants & we'll report a failure
if [ "$VARS_REPORTED" = "0" ]; then
  report_empty_results
  clean_up_early_intermeds "$SUBJECT_FILE" "${CONTROLS[@]}"
fi

# Check for empty lines
echo "Checking vcf formatting..."
bash ${POST_PROCESS_DIR}remove_no_genotype.sh $TEMP_FINAL_VCF "final_no_gx.vcf"
bash ${POST_PROCESS_DIR}remove_no_genotype.sh $TEMP_PREFILTERED_VCF "prefiltered_no_gx.vcf"
rm $TEMP_FINAL_VCF
rm $TEMP_PREFILTERED_VCF
mv "final_no_gx.vcf.gz" $TEMP_FINAL_VCF
mv "prefiltered_no_gx.vcf.gz" $TEMP_PREFILTERED_VCF

# Sort
echo "Sorting..."
$bcftools sort $TEMP_FINAL_VCF | bgzip > "sorted.${TEMP_FINAL_VCF}"
# TODO: when fix formatting on prefiltered vcf, comment two lines below back in
#$bcftools sort $TEMP_PREFILTERED_VCF | bgzip > "sorted.${TEMP_PREFILTERED_VCF}"

rm $TEMP_FINAL_VCF
#rm $TEMP_PREFILTERED_VCF
$bcftools index "sorted.$TEMP_FINAL_VCF"

# Remove coinheriteds
echo "Removing coinheriteds..."
IFS=$','
CONTROL_STRING="${CONTROLS[*]}"
COINHERITED_REMOVED_VCF="coinherited_removed.vcf.gz"
bash ${POST_PROCESS_DIR}remove_coinheriteds.sh "$REFERENCE" "sorted.${TEMP_FINAL_VCF}" "$COINHERITED_REMOVED_VCF" "$SOURCE_DIR" "$CONTROL_STRING"

# Add HD_AF field
echo "Adding kmer-based allele frequencies..." 
AF_ADDED_VCF="hd_af.${COINHERITED_REMOVED_VCF}"
SUBJECT_SAMPLE_NAME=$($bcftools view -h $COINHERITED_REMOVED_VCF | tail -n 1 | awk -F'\t' '{ print $10 }')
bash ${POST_PROCESS_DIR}add_hd_med.add_hd_af.sh "$COINHERITED_REMOVED_VCF" "$SUBJECT_SAMPLE_NAME"
$bcftools index $AF_ADDED_VCF 

# Compose final vcfs
SUBJECT_STRING=$(basename $SUBJECT_FILE)
FINAL_VCF="RUFUS.Final.${SUBJECT_STRING}.combined.vcf"
PREFILTERED_VCF="RUFUS.Prefiltered.${SUBJECT_STRING}.combined.vcf"

# Inject RUFUS command into header
echo "Composing final vcfs..."
$bcftools view -h $AF_ADDED_VCF | head -n -1 > $FINAL_VCF
cat /mnt/rufus.cmd >> $FINAL_VCF
$bcftools view -h $AF_ADDED_VCF | tail -n 1 >> $FINAL_VCF
$bcftools view -H $AF_ADDED_VCF >> $FINAL_VCF
bgzip $FINAL_VCF
$bcftools index "$FINAL_VCF.gz"

#TODO: Comment back in after prefiltered vcf cleaned up
#$bcftools view -h $TEMP_PREFILTERED_VCF | head -n -1 > $PREFILTERED_VCF
#cat /mnt/rufus.cmd >> $PREFILTERED_VCF
#$bcftools view -h $TEMP_PREFILTERED_VCF | tail -n 1 >> $PREFILTERED_VCF
#$bcftools view -H $TEMP_PREFILTERED_VCF >> $PREFILTERED_VCF
#bgzip $PREFILTERED_VCF
#$bcftools index "$PREFILTERED_VCF.gz"
#mv "$PREFILTERED_VCF.gz"* rufus_supplementals/

# Only need to move and rename if did a windowed run
if [ "$WINDOW_SIZE" != "0" ]; then
	mv $TEMP_PREFILTERED_VCF prefiltered.vcf.gz
	mv $TEMP_PREFILTERED_VCF.tbi prefiltered.vcf.gz.tbi
	mv prefiltered.vcf.gz* rufus_supplementals/
fi


# TODO: Separate SVs and SNV/Indels
#echo "Separating snvs/indels and SVs..."

# Cleanup
echo "Cleaning up intermediate post-processing files..."
#rm $TEMP_PREFILTERED_VCF*
rm $TEMP_FINAL_VCF*
#rm "sorted.$TEMP_PREFILTERED_VCF"*
rm "sorted.$TEMP_FINAL_VCF"*
rm $COINHERITED_REMOVED_VCF*
rm "normed.sorted.$TEMP_FINAL_VCF"*
rm -r "/mnt/Intermediates"
rm -r "/mnt/TempOverlap"
rm "/mnt/rufus.cmd"
rm "$AF_ADDED_VCF"*

# Combining supplementals
SUPPLEMENTAL_DIR=/mnt/rufus_supplementals/
# TODO: only do this if not reporting in developer mode
ls ${SUPPLEMENTAL_DIR}*generator.V2.overlap.hashcount.fastq.bam | xargs samtools merge ${SUPPLEMENTAL_DIR}unique_contigs.bam
ls ${SUPPLEMENTAL_DIR}*generator.Mutations.fastq.bam | xargs samtools merge ${SUPPLEMENTAL_DIR}unique_reads.bam
samtools sort ${SUPPLEMENTAL_DIR}unique_contigs.bam -o ${SUPPLEMENTAL_DIR}unique_contigs.sorted.bam
samtools sort ${SUPPLEMENTAL_DIR}unique_reads.bam -o ${SUPPLEMENTAL_DIR}unique_reads.sorted.bam
rm ${SUPPLEMENTAL_DIR}unique_contigs.bam
rm ${SUPPLEMENTAL_DIR}unique_reads.bam
rm ${SUPPLEMENTAL_DIR}*generator.V2.overlap.hashcount.fastq.bam*
rm ${SUPPLEMENTAL_DIR}*generator.Mutations.fastq.bam*

cat ${SUPPLEMENTAL_DIR}*.HashList > ${SUPPLEMENTAL_DIR}unique_kmer_counts.txt
rm ${SUPPLEMENTAL_DIR}*.HashList

clean_up_early_intermeds() "$SUBJECT_FILE" "${CONTROLS[@]}"

echo "Post-processing complete."
end_time=$(date +"%s")
time_delta=$(( $end_time - $start_time ))
hours=$(( time_delta / 3600 ))
minutes=$(( (time_delta % 3600) / 60 ))
seconds=$(( time_delta % 60 ))
printf "RUFUS call stage completed in: %02d:%02d:%02d\n" $hours $minutes $seconds
