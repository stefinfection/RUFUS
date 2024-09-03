#!/bin/bash
#SBATCH --job-name=combine_job     # Job name
#SBATCH --output=combine_job.out   # Output file name
#SBATCH --error=combine_job.err    # Error file name
#SBATCH --time=00:05:00            # Maximum time
#SBATCH --ntasks=1                 # Number of tasks
#SBATCH --cpus-per-task=1          # Number of CPU cores per task
#SBATCH --mem=1G                   # Memory per job

usage() {
	echo "Usage: $0 [-w window_size] [-r reference] -[c control1,control2,control3...] [-h]"
	echo "Options:"
	echo " -w window_size	Required: The size of the window used in the RUFUS run"
	echo " -r reference	Required: The reference used in the RUFUS run"
	echo " -c controls	Required: The control bam files used in the RUFUS run"
	echo " -s subject_string	Required: The name of the subject file"
	echo " -d source_dir	Required: The source directory where the RUFUS vcf(s) are located"
	echo " -h help	Print help message"
	exit 1
}

# initialize vars
controls=()
window_size=0
reference=""
subject_string=""
source_dir=""

# parse command line arguments
while getopts ":w:r:c:s:h" option; do 
	case $option in 
		h) usage;;
		w) window_size=$OPTARG;;
		r) reference=$OPTARG;;
		s) subject_string=$OPTARG;;
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

# trim and combine runs
# need chunk size, path to dirs, name of prefiltered vcf + final vcf
# todo: include header starter file in container

#todo: left off here - get this running locally
bash trim_and_combine.sh $SUBJECT_STRING $

# check for empty lines

# sort

# need final vcf and control bams
# remove inheriteds

# separate SVs and SNV/indels

