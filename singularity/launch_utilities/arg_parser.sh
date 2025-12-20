#!/bin/bash

# Statics
DEFAULT_1MB_CPUS_PER_JOB="12"
DEFAULT_1MB_MEM_PER_JOB="20G"
DEFAULT_WG_CPUS_PER_JOB="40"
DEFAULT_WG_MEM_PER_JOB="150G"


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
  echo "-l slurm_job_array_limit    The maximum amount of jobs slurm allows in an array"
  echo "Optional Arguments:"
  echo "-m kmer_depth_cutoff	The amount of kMers that must overlap the variant to be included in the final call set"
  echo "-w window_size	The size of the windows to run RUFUS on, in units of kilabases (KB); allowed range between 500-5000; defaults to single run of entire genome if not provided"
  echo "-f reference_hash   Jhash file containing reference kMer hash list"
  echo "-x exclude_hash     Single or comma-delimited list of Jhash file(s) containing kMers to exclude from unique hash list"
  echo "-x1kg exclude 1000g population variants     Options are 10 or 100; excludes heterozygous and homozygous variants that are present in 10 or 100 individuals in version 3 of the 1000G project"
  echo "-y path_to_rufus_container   If not provided, will look in current directory for rufus.sif"
  echo "-z rufus_threads	Number of threads provided to RUFUS; defaults to 36 for entire genome; 10 for 1MB windows (NOTE: must be less than cpus_per_call)"
  echo "-e email  The email address to notify with slurm updates"
  echo "-q slurm_job_queue_limit    The maximum amount of jobs able to be ran at once; defaults to 20"
  echo "-t slurm_time_limit   The maximum amount of time to let the slurm job run; defaults to 7 days for full run, or one hour per window (DD-HH:MM:SS)"
  echo "-M memory_per_call    How much memory to allot to the rufus calling stage job; default 150G for entire genome; 20G for 1MB windows (e.g. 150G or 20G)"
  echo "-C cpus_per_call      How many cpus to allot to each rufus calling stage job; default 40 for entire genome; 12 for 1MB windows"
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
SLURM_JOB_LIMIT_RUFUS_ARG="20"
SLURM_ARRAY_JOB_LIMIT_RUFUS_ARG="1000"
SLURM_TIME_LIMIT_RUFUS_ARG=""
CONTAINER_PATH_RUFUS_ARG=""
THREAD_LIMIT_RUFUS_ARG=""
EXCLUDE_HASH_LIST_RUFUS_ARG=()
KG1_EXCLUSION_THRESHOLD="0"
REFERENCE_HASH_RUFUS_ARG=""
MEM_PER_JOB=""
CPUS_PER_JOB=""

# Parse command line options using getopts
while getopts ":d:s:c:b:a:p:r:m:w:e:l:q:t:f:x:y:z:h:M:C" opt; do
    case ${opt} in
        d)
            HOST_DATA_DIR_RUFUS_ARG=$OPTARG
            ;;
        s)
            SUBJECT_RUFUS_ARG=$OPTARG
            ;;
        c)
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
        q)
            SLURM_JOB_LIMIT_RUFUS_ARG=$OPTARG
            ;;
        l)
            SLURM_ARRAY_JOB_LIMIT_RUFUS_ARG=$OPTARG
            ;;
        t)
            SLURM_TIME_LIMIT_RUFUS_ARG=$OPTARG
            ;;
        y)
            CONTAINER_PATH_RUFUS_ARG=$OPTARG
            ;;
		x)
            IFS=',' read -r -a EXCLUDE_HASH_LIST_RUFUS_ARG <<< "$OPTARG"
			;;
        x1kg)
            KG1_EXCLUSION_THRESHOLD=$OPTARG
            ;;
		f)
			REFERENCE_HASH_RUFUS_ARG=$OPTARG
		    ;;
		z)
			THREAD_LIMIT_RUFUS_ARG=$OPTARG
			;;
        M)
			MEM_PER_JOB=$OPTARG
			;;
        C)
			CPUS_PER_JOB=$OPTARG
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
if [[ -z "$HOST_DATA_DIR_RUFUS_ARG" || -z "$SUBJECT_RUFUS_ARG" || -z "$GENOME_BUILD_RUFUS_ARG" || -z "$REFERENCE_RUFUS_ARG" || -z "$SLURM_ACCOUNT_RUFUS_ARG" || -z "$SLURM_PARTITION_RUFUS_ARG" || -z "$SLURM_ARRAY_JOB_LIMIT_RUFUS_ARG" ]]; then
    echo "ERROR: Missing required argument(s); please see usage instructions with -h." >&2
fi

# Add on a trailing slash to dir, just in case user omits
HOST_DATA_DIR_RUFUS_ARG="${HOST_DATA_DIR_RUFUS_ARG}/"

# Check that data directory exists
if [ ! -d "$HOST_DATA_DIR_RUFUS_ARG" ]; then
	echo "ERROR: provided data_directory argument is not a directory." >&2
fi

# Check that subject file is in provided data directory
if [ ! -f "${HOST_DATA_DIR_RUFUS_ARG}${SUBJECT_RUFUS_ARG}" ]; then
	echo "ERROR: provided subject file $SUBJECT_RUFUS_ARG does not exist in the provided data directory or cannot be read." >&2
fi

# Check that all of the control files are in the provided data directory
for control in "${CONTROLS_RUFUS_ARG[@]}"; do
	if [ ! -f "${HOST_DATA_DIR_RUFUS_ARG}${control}" ]; then
		echo "ERROR: provided control file $control does not exist in the provided data directory or cannot be read." >&2
	else
		if [ "$CONTROL_STRING_RUFUS_ARG" == "" ]; then
			CONTROL_STRING_RUFUS_ARG="$control"
		else
			CONTROL_STRING_RUFUS_ARG="${CONTROL_STRING_RUFUS_ARG}, $control"
		fi
	fi
done

# Check that reference file is in provided data directory
if [ ! -f "${HOST_DATA_DIR_RUFUS_ARG}${REFERENCE_RUFUS_ARG}" ]; then
	echo "ERROR: provided reference file $REFERENCE_RUFUS_ARG does not exist in the provided data directory or cannot be read." >&2
    exit 1
fi

# Check that window size is in valid range
# Check if time limit has been assigned, if not - use defaults for full mode or windowed mode
if [ "$WINDOW_SIZE_RUFUS_ARG" -eq 0 ]; then
	if [ -z "$SLURM_TIME_LIMIT_RUFUS_ARG" ]; then
		SLURM_TIME_LIMIT_RUFUS_ARG="7-00:00:00"
	fi
    CPUS_PER_JOB=${CPUS_PER_JOB:-$DEFAULT_WG_CPUS_PER_JOB}
    MEM_PER_JOB=${MEM_PER_JOB:-$DEFAULT_WG_MEM_PER_JOB}
    THREAD_LIMIT_RUFUS_ARG=${THREAD_LIMIT_RUFUS_ARG:-36}
elif [ "$WINDOW_SIZE_RUFUS_ARG" -ne 1000 ]; then
	echo "ERROR: only windows of 1000 (1MB) supported currently" >&2
    exit 1
else
	if [ -z $SLURM_TIME_LIMIT_RUFUS ]; then
		SLURM_TIME_LIMIT_RUFUS_ARG="01:00:00"
	fi
    CPUS_PER_JOB=${CPUS_PER_JOB:-$DEFAULT_1MB_CPUS_PER_JOB}
    MEM_PER_JOB=${MEM_PER_JOB:-$DEFAULT_1MB_MEM_PER_JOB}
    THREAD_LIMIT_RUFUS_ARG=${THREAD_LIMIT_RUFUS_ARG:-10}
fi

echo "debugging: cpus per job: $CPUS_PER_JOB and thread limit $THREAD_LIMIT_RUFUS_ARG" >&2

if [ "$THREAD_LIMIT_RUFUS_ARG" -ge "$CPUS_PER_JOB" ]; then
	echo "ERROR: thread limit ($THREAD_LIMIT_RUFUS_ARG) must be less than cpus per job ($CPUS_PER_JOB)." >&2
	exit 1
fi
# TODO: have to check that window size is either 0 or 1000 because no other mini-hashes supported yet
# Check that 1000kg exclusion threshold is valid
if [ "$KG1_EXCLUSION_THRESHOLD" -ne 10 ] && [ "$KG1_EXCLUSION_THRESHOLD" -ne 100 ] && [ "$KG1_EXCLUSION_THRESHOLD" -ne 0 ]; then
    echo "ERROR: 1000kg exclusion threshold must be either 10 or 100"
elif [ "$KG1_EXCLUSION_THRESHOLD" -ne 0 ] && ([ "$WINDOW_SIZE_RUFUS_ARG" -ne 1000 ]); then
    echo "ERROR: 1000kg exclusion hashes can only be used with full genome run or 1000kb window size"
fi

# Check that if path to image not provided, it's in the current dir
if [ -z $CONTAINER_PATH_RUFUS_ARG ]; then
	if [ ! -f "rufus.sif" ]; then
		echo "ERROR: rufus.sif not in current directory - please provide path to container or put it in this one under rufus.sif"
	fi
fi

# Export variables for use in the main script
export HOST_DATA_DIR_RUFUS_ARG
export SUBJECT_RUFUS_ARG
export CONTROL_STRING_RUFUS_ARG
export CONTROLS_RUFUS_ARG
export GENOME_BUILD_RUFUS_ARG
export SLURM_ACCOUNT_RUFUS_ARG
export SLURM_PARTITION_RUFUS_ARG
export REFERENCE_RUFUS_ARG
export KMER_DEPTH_CUTOFF_RUFUS_ARG
export WINDOW_SIZE_RUFUS_ARG
export EMAIL_RUFUS_ARG
export SLURM_JOB_LIMIT_RUFUS_ARG
export SLURM_ARRAY_JOB_LIMIT_RUFUS_ARG
export SLURM_TIME_LIMIT_RUFUS_ARG
export CONTAINER_PATH_RUFUS_ARG
export THREAD_LIMIT_RUFUS_ARG
export EXCLUDE_HASH_LIST_RUFUS_ARG
export KG1_EXCLUSION_THRESHOLD
export REFERENCE_HASH_RUFUS_ARG
export MEM_PER_JOB
export CPUS_PER_JOB