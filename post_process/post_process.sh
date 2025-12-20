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

report_empty_and_exit() {
  echo "RUFUS did not find any variants for the provided parameters. Please adjust and try again."
  echo "RUFUS did not find any variants for the provided parameters. Please adjust and try again." >&2
  echo "RUFUS did not find any variants for the provided parameters. Please adjust and try again." > results.out
  exit 0
}

# Cleans up intermediate files, reports no variants found in both out + error, and exits failure code
clean_up_post_temps() {
  files=("$TEMP_FINAL_VCF" "sorted.$TEMP_FINAL_VCF" \
  "$COINHERITED_REMOVED_VCF" "normed.sorted.$TEMP_FINAL_VCF" \
  "$AF_ADDED_VCF" "rufus.cmd" "final_no_gx.vcf" )

  for file in "${files[@]}"; do
    find . -type f -name "$file" -delete
  done

  find /mnt -type d -name "Intermediates" -delete 
  find /mnt -type d -name "TempOverlap" -delete 
}
trap 'clean_up_post_temps' EXIT

# We don't want to do this unless post-processing completes without error
clean_up_calls() {
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

if [[ -s "$REFERENCE" ]]; then
	    echo "ERROR: Missing required option -r (reference)" >&2
fi

if [[ -d "$SOURCE_DIR" ]]; then
	echo "ERROR: Missing required option -d (source directory for RUFUS vcf(s))" >&2
fi

if [[ -s "$SUBJECT_FILE" ]]; then
	echo "ERROR: Missing required option -s (subject cram/bam)" >&2
fi

if [ ${#CONTROLS[@]} -eq 0 ]; then
	    echo "ERROR: Must supply at least one control cram/bam" >&2
fi


cd "$SOURCE_DIR" || exit 1
echo "RUFUS post-process version E-0.1.0"
date
start_time=$(date +"%s")

POST_PROCESS_DIR=/opt/RUFUS/post_process/
TEMP_FINAL_VCF="temp.RUFUS.Final.${SUBJECT_FILE}.combined.vcf.gz"
TEMP_PREFILTERED_VCF="temp.RUFUS.Prefiltered.${SUBJECT_FILE}.combined.vcf.gz"
GERMLINE_VCF="with_germline.RUFUS.Final.${SUBJECT_FILE}.combined.vcf.gz"

# Slight name change if not doing a windowed run
if [ "$WINDOW_SIZE" = "0" ]; then
	TEMP_FINAL_VCF="temp.RUFUS.Final.${SUBJECT_FILE}.vcf.gz"
fi

# Check to see if temp vcf(s) exists, if not report empty results and exit
if read -r first_match < <(compgen -G "temp.RUFUS.Final*vcf.gz"); then
    echo "Found temporary vcf(s)"
else
  report_empty_and_exit
fi

# If windowed mode, trim and combine region
if [ "$WINDOW_SIZE" != "0" ]; then
	IFS=$'\t'
	echo "Windowed run performed, trimming and combining region vcfs..."
	bash ${POST_PROCESS_DIR}trim_and_combine.sh "$SUBJECT_FILE"  "$WINDOW_SIZE" "${CONTROLS[@]}"
fi

# Get number of variants reported
VARS_REPORTED=$($bcftools view -H "$TEMP_FINAL_VCF" | wc -l)

# Check for empty vcf AFTER trimming and combining
# If we don't have any variants here, the entire run didn't find any variants & we'll report a failure
if [ "$VARS_REPORTED" = "0" ]; then
  clean_up_early_intermeds "$SUBJECT_FILE" "${CONTROLS[@]}"
  report_empty_and_exit
fi

# Check for empty lines
echo "Checking vcf formatting..."
bash ${POST_PROCESS_DIR}remove_no_genotype.sh "$TEMP_FINAL_VCF" "final_no_gx.vcf"
mv "final_no_gx.vcf.gz" "$TEMP_FINAL_VCF"

# Sort
echo "Sorting..."
$bcftools sort "$TEMP_FINAL_VCF" | bgzip > "sorted.${TEMP_FINAL_VCF}"
$bcftools index "sorted.$TEMP_FINAL_VCF"

# Remove coinheriteds
echo "Removing coinheriteds..."
IFS=$','
COINHERITED_REMOVED_VCF="coinherited_removed.vcf.gz"
bash ${POST_PROCESS_DIR}remove_coinheriteds.sh "$REFERENCE" "sorted.${TEMP_FINAL_VCF}" "$COINHERITED_REMOVED_VCF" "$SOURCE_DIR" "$WINDOW_SIZE" "${CONTROLS[@]}"

# Add HD_AF field
echo "Adding kmer-based allele frequencies..." 
AF_ADDED_VCF="hd_af.${COINHERITED_REMOVED_VCF}"
SUBJECT_SAMPLE_NAME=$($bcftools view -h $COINHERITED_REMOVED_VCF | tail -n 1 | awk -F'\t' '{ print $10 }')
bash ${POST_PROCESS_DIR}add_hd_med.add_hd_af.sh "$COINHERITED_REMOVED_VCF" "$SUBJECT_SAMPLE_NAME"
$bcftools index $AF_ADDED_VCF 

# Compose final vcfs
SUBJECT_STRING=$(basename "$SUBJECT_FILE")
FINAL_VCF="RUFUS.Final.${SUBJECT_STRING}.combined.vcf"
#PREFILTERED_VCF="RUFUS.Prefiltered.${SUBJECT_STRING}.combined.vcf"

# Inject RUFUS command into header
echo "Composing final vcfs..."
$bcftools view -h $AF_ADDED_VCF | head -n -1 > "$FINAL_VCF"
cat /mnt/rufus.cmd >> "$FINAL_VCF"
$bcftools view -h $AF_ADDED_VCF | tail -n 1 >> "$FINAL_VCF"
$bcftools view -H $AF_ADDED_VCF >> "$FINAL_VCF"
bgzip "$FINAL_VCF"
$bcftools index "$FINAL_VCF.gz"

#TODO: Comment back in after prefiltered vcf cleaned up
# $bcftools view -h $TEMP_PREFILTERED_VCF | head -n -1 > $PREFILTERED_VCF
# cat /mnt/rufus.cmd >> $PREFILTERED_VCF
# $bcftools view -h $TEMP_PREFILTERED_VCF | tail -n 1 >> $PREFILTERED_VCF
# $bcftools view -H $TEMP_PREFILTERED_VCF >> $PREFILTERED_VCF
# bgzip $PREFILTERED_VCF
# $bcftools index "$PREFILTERED_VCF.gz"
# mv "$PREFILTERED_VCF.gz"* rufus_supplementals/

# Only need to move and rename if did a windowed run
# if [ "$WINDOW_SIZE" != "0" ]; then
# 	mv $TEMP_PREFILTERED_VCF prefiltered.vcf.gz
# 	mv $TEMP_PREFILTERED_VCF.tbi prefiltered.vcf.gz.tbi
# 	mv prefiltered.vcf.gz* rufus_supplementals/
# fi


# TODO: Separate SVs and SNV/Indels
#echo "Separating snvs/indels and SVs..."

# Combining supplementals
#SUPPLEMENTAL_DIR=/mnt/rufus_supplementals/
# TODO: only do this if not reporting in developer mode
# ls ${SUPPLEMENTAL_DIR}*generator.V2.overlap.hashcount.fastq.bam | xargs samtools merge ${SUPPLEMENTAL_DIR}unique_contigs.bam
# ls ${SUPPLEMENTAL_DIR}*generator.Mutations.fastq.bam | xargs samtools merge ${SUPPLEMENTAL_DIR}unique_reads.bam
# samtools sort ${SUPPLEMENTAL_DIR}unique_contigs.bam -o ${SUPPLEMENTAL_DIR}unique_contigs.sorted.bam
# samtools sort ${SUPPLEMENTAL_DIR}unique_reads.bam -o ${SUPPLEMENTAL_DIR}unique_reads.sorted.bam
# rm ${SUPPLEMENTAL_DIR}unique_contigs.bam
# rm ${SUPPLEMENTAL_DIR}unique_reads.bam
# rm ${SUPPLEMENTAL_DIR}*generator.V2.overlap.hashcount.fastq.bam*
# rm ${SUPPLEMENTAL_DIR}*generator.Mutations.fastq.bam*

# cat ${SUPPLEMENTAL_DIR}*.HashList > ${SUPPLEMENTAL_DIR}unique_kmer_counts.txt
# rm ${SUPPLEMENTAL_DIR}*.HashList

clean_up_calls

echo "Post-processing complete."
end_time=$(date +"%s")
time_delta=$(( end_time - start_time ))
hours=$(( time_delta / 3600 ))
minutes=$(( (time_delta % 3600) / 60 ))
seconds=$(( time_delta % 60 ))
printf "RUFUS call stage completed in: %02d:%02d:%02d\n" $hours $minutes $seconds