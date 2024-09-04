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
controls=()
window_size=0
reference=""
subject_file=""
source_dir=""

# parse command line arguments
while getopts ":w:r:c:s:d:h" option; do 
	case $option in 
		h) usage;;
		w) window_size=$OPTARG;;
		r) reference=$OPTARG;;
		s) subject_file=$OPTARG;;
		d) source_dir=$OPTARG;;
		c) IFS=',' read -r -a controls <<< "$OPTARG";;
		\?) echo "Invalid option: -$OPTARG" >&2
		    usage;;	
		*) echo "Option -$OPTARG requires an argument" >&2
		   usage;;
	esac
done
shift $((OPTIND-1))

# check for mandatory command line arguments
if [[ -z "$window_size" ]]; then
	    echo "Error: Missing required option -w (window size)" >&2
	        usage
fi

if [[ -z "$reference" ]]; then
	    echo "Error: Missing required option -r (reference)" >&2
	        usage
fi

if [[ -z "$source_dir" ]]; then
	echo "Error: Missing required option -d (source directory for RUFUS vcf(s))" >&2
	        usage
fi

if [ ${#controls[@]} -eq 0 ]; then
	    echo "Error: Must supply at least one control bam" >&2
	        usage
fi

echo "RUFUS post-process version C.0.1"
echo "RUFUS post-process command was: postProcess.sh $@"
date
echo "RUFUS_$0 $@" >> .rufus.cmd

cd $source_dir

subject_string=$(basename $subject_file)
POST_PROCESS_DIR=/opt/RUFUS/post_process/

TEMP_FINAL_VCF="temp.RUFUS.Final.${subject_string}.combined.vcf.gz"
TEMP_PREFILTERED_VCF="temp.RUFUS.Prefiltered.${subject_string}.combined.vcf.gz"

if [ "$window_size" != "0" ]; then
	echo "Windowed run performed, trimming and combining region vcfs..."
	IFS=$'\t'
	control_string="${controls[*]}"
	bash ${POST_PROCESS_DIR}trim_and_combine.sh $subject_string $control_string $source_dir $window_size
	mv $PREFILTERED_VCF rufus_supplementals/
fi

# Check for empty lines
bash ${POST_PROCESS_DIR}remove_no_genotype.sh $TEMP_FINAL_VCF
bash ${POST_PROCESS_DIR}remove_no_genotype.sh $TEMP_PREFILTERED_VCF

# Sort
bcftools sort $TEMP_FINAL_VCF > "sorted.${TEMP_FINAL_VCF}"
bcftools sort $TEMP_PREFILTERED_VCF > "sorted.${TEMP_PREFILTERED_VCF}"

# Remove coinheriteds
# TODO: left off here - we are already inside of a batch script, don't need to make another one - call just a .sh here instead


# Separate SVs and SNV/Indels


FINAL_VCF="RUFUS.Final.${subject_string}.combined.vcf.gz"
PREFILTERED_VCF="RUFUS.Prefiltered.${subject_string}.combined.vcf.gz"

# Inject RUFUS command into header
bcftools view -h $TEMP_FINAL_VCF | head -n -1 > $FINAL_VCF
cat .rufus.cmd >> $FINAL_VCF
bcftools view -h $TEMP_FINAL_VCF | tail -n 1 >> $FINAL_VCF
bcftools view -H $TEMP_FINAL_VCF >> $FINAL_VCF

bcftools view -h $TEMP_PREFILTERED_VCF | head -n -1 > $PREFILTERED_VCF
cat .rufus.cmd >> $PREFILTERED_VCF
bcftools view -h $TEMP_PREFILTERED_VCF | tail -n 1 >> $PREFILTERED_VCF
bcftools view -H $TEMP_PREFILTERED_VCF >> $PREFILTERED_VCF

# Cleanup
mv $PREFILTERED_VCF rufus_supplementals/
rm $TEMP_PREFILTERED_VCF
rm $TEMP_FINAL_VCF
rm "sorted.$TEMP_PREFILTERED_VCF"
rm "sorted.$TEMP_FINAL_VCF"

