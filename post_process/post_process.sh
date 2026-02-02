#!/bin/bash

usage() {
	echo "Usage: $0 [-s subject]-d source_dir] [-h]"
	echo "Options:"
	echo " -s subject_file	Required: The name of the subject file: must be the same as that supplied to the RUFUS run"
	echo " -d source_dir	Required: The source directory where the RUFUS vcf(s) are located" # TODO: make this not required
	echo " -h help	Print help message"
	exit 1
}

# Cleans up intermediate files, reports no variants found in both out + error, and exits failure code
clean_up_post_temps() {
  files=("$TEMP_FINAL_VCF" "rufus.cmd" )

  for file in "${files[@]}"; do
	if [ "$file" != "" ]; then
    	find . -maxdepth 1 -type f -name "$file*" -delete
	fi
  done

  find . -type d -name "Intermediates" -exec rm -rf {} +
  find . -type d -name "TempOverlap" -exec rm -rf {} +
}
trap 'clean_up_post_temps' EXIT

# parse command line arguments
SUBJECT_FILE=""
SOURCE_DIR="."

while getopts "h:s:d:" option; do 
	case $option in 
		h) usage;;
		s) SUBJECT_FILE=$OPTARG;;
		d) SOURCE_DIR=$OPTARG;;
		\?) echo "Invalid option: -$OPTARG" >&2
		    usage;;	
		*) echo "Option -$OPTARG requires an argument" >&2
		   usage;;
	esac
done
shift $((OPTIND-1))

if [[ -z "$SUBJECT_FILE" ]]; then
	echo "ERROR: Missing required option -s (subject cram/bam)" >&2
	exit 1
fi


echo "RUFUS post-process version E-0.1.0"
date
start_time=$(date +"%s")

# Have to define all of these before first possible exit
SUBJECT_STRING=$(basename "$SUBJECT_FILE")
TEMP_FINAL_VCF="temp.RUFUS.Final.${SUBJECT_STRING}.combined.vcf.gz"
TEMP_PREFILTERED_VCF="temp.RUFUS.Prefiltered.${SUBJECT_STRING}.combined.vcf.gz"
FINAL_VCF="RUFUS.Final.${SUBJECT_STRING}.vcf"

# Slight name change if not doing a windowed run
if [ "$WINDOW_SIZE" -eq 0 ]; then
	TEMP_FINAL_VCF="temp.RUFUS.Final.${SUBJECT_FILE}.vcf.gz"
fi

# Check to see if temp vcf(s) exists, if not report empty results and exit
if read -r first_match < <(compgen -G "temp.RUFUS.Final*vcf.gz"); then
    echo "Found temporary vcf(s)"
else
  	echo "Could not find any temporary vcf(s) from calling stage. Exiting..."
fi

MERGED="merged.vcf"
NO_HEAD="no_header.vcf"
FINAL_GZ="${FINAL_VCF}.gz"

# Concat vcfs if in windowed mode
if [ "$WINDOW_SIZE" -ne 0 ]; then
	# Concatenate all intermediates
	find . -maxdepth 1 -type f -name 'temp.RUFUS.Final.*.vcf.gz' -print |
	sort -V > regions.txt

	split -l 200 regions.txt regions.chunk.

	for f in regions.chunk.*; do
			echo "-- BCFTools ---------------------------"
			bcftools concat -a -D -f $f -Oz -o "$f.vcf.gz"
			bcftools index "$f.vcf.gz"
	done

	ls regions.chunk.*.vcf.gz > final.list
	bcftools concat -a -D -f final.list -Ov -o $MERGED

	bcftools sort -T "tmp_bcftools.XXXXXX" -O -o "$NO_HEAD" "$MERGED" \
	|| { echo "Error: bcftools sort failed"; exit 1; }

else
	# Just sort and add header for whole genome mode
	bcftools sort -T "tmp_bcftools.XXXXXX" -O -o "$NO_HEAD" "$TEMP_FINAL_VCF" \
		|| { echo "Error: bcftools sort failed"; exit 1; }
fi

# Inject rufus command into header
bcftools view -h $NO_HEAD | head -n -1 > "$FINAL_VCF" || { echo "Error: bcftools view on merged vcf failed"; exit 1; }
cat rufus.cmd >> "$FINAL_VCF" || { echo "Could not find rufus.cmd"; exit 1; }
bcftools view -h $NO_HEAD | tail -n 1 >> "$FINAL_VCF"
bcftools view -H $NO_HEAD >> "$FINAL_VCF" 

bgzip "$FINAL_VCF" || { echo "Error: bgzip failed on final vcf"; exit 1; }
tabix -f -p vcf "$FINAL_GZ" \
|| { echo "Error: tabix failed"; exit 1; }
echo "-----------------------------------------"
echo "Done: $FINAL_GZ"

echo "Concatenating & sorting complete."
end_time=$(date +"%s")
time_delta=$(( end_time - start_time ))
hours=$(( time_delta / 3600 ))
minutes=$(( (time_delta % 3600) / 60 ))
seconds=$(( time_delta % 60 ))
printf "RUFUS combine stage completed in: %02d:%02d:%02d\n" $hours $minutes $seconds