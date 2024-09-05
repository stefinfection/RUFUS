#!/bin/bash

usage() {
	echo "Usage: $0 [-w window_size] [-r reference] -[c control1,control2,control3...] [-h]"
	echo "Options:"
	echo " -w window_size	Required: The size of the window used in the RUFUS run"
	echo " -r reference	Required: The reference used in the RUFUS run"
	echo " -c controls	Required: The control bam files used in the RUFUS run"
	echo " -s subject_file	Required: The name of the subject file: must be the same as that supplied to the RUFUS run"
	echo " -d source_dir	Required: The source directory where the RUFUS vcf(s) are located"
	echo " -h help	Print help message"
	exit 1
}

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
	    echo "Error: Missing required option -w (window size)" >&2
	        usage
fi

if [[ -z "$REFERENCE" ]]; then
	    echo "Error: Missing required option -r (reference)" >&2
	        usage
fi

if [[ -z "$SOURCE_DIR" ]]; then
	echo "Error: Missing required option -d (source directory for RUFUS vcf(s))" >&2
	        usage
fi

if [ ${#CONTROLS[@]} -eq 0 ]; then
	    echo "Error: Must supply at least one control bam" >&2
	        usage
fi

cd $SOURCE_DIR
echo "RUFUS post-process version C.0.1"
date
echo "RUFUS_$0 $@" >> rufus.cmd

SUBJECT_STRING=$(basename $SUBJECT_FILE)
POST_PROCESS_DIR=/opt/RUFUS/post_process/

TEMP_FINAL_VCF="temp.RUFUS.Final.${SUBJECT_STRING}.combined.vcf.gz"
TEMP_PREFILTERED_VCF="temp.RUFUS.Prefiltered.${SUBJECT_STRING}.combined.vcf.gz"

if [ "$WINDOW_SIZE" != "0" ]; then
	IFS=$'\t'
	TAB_DELIM_CONTROL_STRING="${CONTROLS[*]}"
	echo "Windowed run performed, trimming and combining region vcfs..."
	bash ${POST_PROCESS_DIR}trim_and_combine.sh $SUBJECT_STRING $TAB_DELIM_CONTROL_STRING $WINDOW_SIZE
	mv $PREFILTERED_VCF rufus_supplementals/ # todo: am I doing this in trim and combine also?
fi

# Check for empty lines
echo "Checking vcf formatting..."
bash ${POST_PROCESS_DIR}remove_no_genotype.sh $TEMP_FINAL_VCF
bash ${POST_PROCESS_DIR}remove_no_genotype.sh $TEMP_PREFILTERED_VCF

# Sort
echo "Sorting..."
bcftools sort $TEMP_FINAL_VCF > "sorted.${TEMP_FINAL_VCF}"
bcftools sort $TEMP_PREFILTERED_VCF > "sorted.${TEMP_PREFILTERED_VCF}"

FINAL_VCF="RUFUS.Final.${subject_string}.combined.vcf.gz"
PREFILTERED_VCF="RUFUS.Prefiltered.${subject_string}.combined.vcf.gz"

# Remove coinheriteds
echo "Removing coinheriteds..."
IFS=$','
CONTROL_STRING="${CONTROLS_RUFUS_ARG[*]}"
bash ${POST_PROCESS_DIR}remove_coinheriteds.sh "$REFERENCE" "sorted.${TEMP_FINAL_VCF}" "$FINAL_VCF" "$SOURCE_DIR" "$CONTROL_STRING"

# Inject RUFUS command into header
echo "Composing final vcfs..."
bcftools view -h $COINHERITED_REMOVED_VCF | head -n -1 > $FINAL_VCF
cat rufus.cmd >> $FINAL_VCF
bcftools view -h $COINHERITED_REMOVED_VCF | tail -n 1 >> $FINAL_VCF
bcftools view -H $COINHERITED_REMOVED_VCF >> $FINAL_VCF

bcftools view -h $TEMP_PREFILTERED_VCF | head -n -1 > $PREFILTERED_VCF
cat rufus.cmd >> $PREFILTERED_VCF
bcftools view -h $TEMP_PREFILTERED_VCF | tail -n 1 >> $PREFILTERED_VCF
bcftools view -H $TEMP_PREFILTERED_VCF >> $PREFILTERED_VCF

# TODO: Separate SVs and SNV/Indels
#echo "Separating snvs/indels and SVs..."

# Cleanup
echo "Cleaning up intermediate post-processing files..."
#TODO: cleanup combined vcf
mv $PREFILTERED_VCF rufus_supplementals/
rm $TEMP_PREFILTERED_VCF
rm $TEMP_FINAL_VCF
rm "sorted.$TEMP_PREFILTERED_VCF"
rm "sorted.$TEMP_FINAL_VCF"

echo "Post-processing complete."

