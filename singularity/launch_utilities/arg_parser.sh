#!/bin/bash

usage() {
    echo "Usage: $0 [-s subject] [-c control1,control2,control3...] [-b genome_build] [-a slurm_account] [-p slurm_partition] ...options"
    echo "Required Arguments:"
	echo "-d data_directory	The directory containing the subject, control, and reference files to be used in the run"
    echo "-s subject    The subject sample of interest; must be located in data_directory"
    echo "-c control(s) A single control or comma-delimited array of multiple controls; must be located in data_directory"
    echo "-b genome_build  The desired genome build; currently only supports GRCh38"
    echo "-r reference  The reference file matching the genome build; must be located in data_directory"
    echo "-a slurm_account  The account for the slurm job"
    echo "-p slurm_partition    The partition for the slurm job"
    echo "Optional Arguments:"
    echo "-r reference  If not provided, will automatically provide GRCh38 file"
	echo "-m kmer_depth_cutoff	The amount of kMers that must overlap the variant to be included in the final call set"
	echo "-w window_size	The size of the windows to run RUFUS on, in units of kilabases (KB); allowed range between 500-5000; defaults to single run of entire genome if not provided" 
    echo "-e email  The email address to notify with slurm updates"
    echo "-l slurm_job_limit    The maximum amount of jobs able to be queued at once"
    echo "-t slurm_time_limit   The maximum amount of time to let the slurm job run; defaults to 7 days (DD-HH:MM:SS)"
    echo "-p path_to_rufus_container    If not provided, will look in current directory for rufus.sif"
	echo "-h help	Print usage"
	exit 1
}

# Initialize variables
HOST_DATA_DIR_RUFUS_ARG=""
SUBJECT_RUFUS_ARG=""
CONTROL_STRING_RUFUS_ARG=""
CONTROLS_RUFUS_ARG=()
GENOME_BUILD_RUFUS_ARG="GRCh38"
SLURM_ACCOUNT_RUFUS_ARG=""
SLURM_PARTITION_RUFUS_ARG=""
REFERENCE_RUFUS_ARG=""
KMER_DEPTH_CUTOFF_RUFUS_ARG="5"
WINDOW_SIZE_RUFUS_ARG="0"
EMAIL_RUFUS_ARG=""
SLURM_JOB_LIMIT_RUFUS_ARG="7-00:00:00"
SLURM_TIME_LIMIT_RUFUS_ARG=""
CONTAINER_PATH_RUFUS_ARG=""
REF_HASH_RUFUS_ARG=""

# Parse command line options using getopts
while getopts ":d:s:c:b:a:p:r:m:w:e:l:t:f:h" opt; do
    case ${opt} in
		d)	
			HOST_DATA_DIR_RUFUS_ARG=$OPTARG
			;;
        s)
            SUBJECT_RUFUS_ARG=$OPTARG
            ;;
        c)
			# todo: need to test this, don't know if I can do two ops on same OPTARG
			CONTROL_STRING_RUFUS_ARG=$OPTARG
            IFS=',' read -r -a CONTROLS_RUFUS_ARG <<< "$OPTARG"
            ;;
        b)
            GENOME_BUILD_RUFUS_ARG=$OPTARG
            ;;
        a)
            SLURM_ACCOUNT_RUFUS_ARG=$OPTARG
            ;;
        p)
            SLURM_PARTITION_RUFUS_ARG=$OPTARG
            ;;
        r)
            REFERENCE_RUFUS_ARG=$OPTARG
            ;;
        m)
            KMER_DEPTH_CUTOFF_RUFUS_ARG=$OPTARG
            ;;
		w)
			WINDOW_SIZE_RUFUS_ARG=$OPTARG
			;;
        e)
            EMAIL_RUFUS_ARG=$OPTARG
            ;;
        l)
            SLURM_JOB_LIMIT_RUFUS_ARG=$OPTARG
            ;;
        t)
            SLURM_TIME_LIMIT_RUFUS_ARG=$OPTARG
            ;;
        f)
            CONTAINER_PATH_RUFUS_ARG=$OPTARG
            ;;
        h)
            usage
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            usage
            ;;
    esac
done
shift $((OPTIND - 1))

# Check for required strings
if [[ -z "$HOST_DATA_DIR_RUFUS_ARG" || -z "$SUBJECT_RUFUS_ARG" || -z "$GENOME_BUILD_RUFUS_ARG" || -z "$REFERENCE_RUFUS_ARG" || -z "$SLURM_ACCOUNT_RUFUS_ARG" || -z "$SLURM_PARTITION_RUFUS_ARG" ]]; then
    echo "Error: Missing required argument(s)." >&2
    usage
fi

# Add on a trailing slash to dir, just in case user omits
HOST_DATA_DIR_RUFUS_ARG="${HOST_DATA_DIR_RUFUS_ARG}/"

# Check that data directory exists
if [ ! -d "$HOST_DATA_DIR_RUFUS_ARG" ]; then
	echo "Error: provided data_directory argument is not a directory." >&2
	usage	
fi

# Check that subject file is in provided data directory
if [ ! -f "${HOST_DATA_DIR_RUFUS_ARG}${SUBJECT_RUFUS_ARG}" ]; then
	echo "Error: provided subject file $SUBJECT_RUFUS_ARG does not exist in the provided data directory or cannot be read." >&2
	usage	
fi

# Check that all of the control files are in the provided data directory
for control in "${CONTROLS_RUFUS_ARG[@]}"; do
	if [ ! -f "${HOST_DATA_DIR_RUFUS_ARG}${control}" ]; then
		echo "Error: provided control file $controls does not exist in the provided data directory or cannot be read." >&2
		usage	
	fi
done

# Check that reference file is in provided data directory
if [ ! -f "${HOST_DATA_DIR_RUFUS_ARG}${REFERENCE_RUFUS_ARG}" ]; then
	echo "Error: provided reference file $REFERENCE_RUFUS_ARG does not exist in the provided data directory or cannot be read." >&2
	usage	
fie

# If build is GRCh38, include prebuilt hash
if [ "$GENOME_BUILD_RUFUS_ARG" = "GRCh38" ]; then
	REF_HASH_RUFUS_ARG="/opt/RUFUS/resources/references/prebuilt_hashes/GRCh38_full_analysis_set_plus_decoy_hla.25.Jhash"
fi

# Export variables for use in the main script
export HOST_DATA_DIR_RUFUS_ARG
export SUBJECT_RUFUS_ARG
export CONTROL_STRING_RUFUS_ARG
export CONTROLS_RUFUS_ARG #TODO: get rid of this - just check here and pass as string
export GENOME_BUILD_RUFUS_ARG
export SLURM_ACCOUNT_RUFUS_ARG
export SLURM_PARTITION_RUFUS_ARG
export REFERENCE_RUFUS_ARG
export KMER_DEPTH_CUTOFF_RUFUS_ARG
export WINDOW_SIZE_RUFUS_ARG
export EMAIL_RUFUS_ARG
export SLURM_JOB_LIMIT_RUFUS_ARG
export SLURM_TIME_LIMIT_RUFUS_ARG
export CONTAINER_PATH_RUFUS_ARG
export REF_HASH_RUFUS_ARG
