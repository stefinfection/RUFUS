#!/bin/bash

# Statics
DEFAULT_1MB_CPUS_PER_JOB="12"
DEFAULT_1MB_MEM_PER_JOB="20G"
DEFAULT_WG_CPUS_PER_JOB="40"
DEFAULT_WG_MEM_PER_JOB="150G"


usage() {
  echo "Usage: $0 [-s subject] [-c control1,control2,control3...] [-b genome_build] [-a slurm_account] [-p slurm_partition] ...options"
  echo "Required Arguments:"
  echo "-s subject    Full path to the subject sample BAM/CRAM"
  echo "-c control(s) A single control or comma-delimited array of multiple controls (full paths)"
  echo "-b genome_build  The desired genome build; currently only supports GRCh38"
  echo "-r reference  Full path to the reference file matching the genome build"
  echo "-a slurm_account  The account for the slurm job"
  echo "-p slurm_partition    The partition for the slurm job"
  echo "-l slurm_job_array_limit    The maximum amount of jobs slurm allows in an array"
  echo "Optional Arguments:"
  echo "-m kmer_depth_cutoff	The amount of kMers that must overlap the variant to be included in the final call set"
  echo "-w window_size	The size of the windows to run RUFUS on, in units of kilabases (KB); allowed range between 500-5000; defaults to single run of entire genome if not provided"
  echo "-f reference_hash   Full path to Jhash file containing reference kMer hash list"
  echo "-x exclude_hash     Single or comma-delimited list of full paths to Jhash file(s) containing kMers to exclude (static, same for all regions)"
  echo "-K kg1_hash_dir     Full path to directory of per-region KG1 Jhash files (files named *{region}*.Jhash)"
  echo "-G kg1_version      KG1 hash version to download from S3 (e.g., v3.0)"
  echo "-D ctrl_hash_dir    Full path to directory of per-region control Jhash files (files named *{region}*.Jhash)"
  echo "-V ctrl_version     Control hash version to download from S3 (e.g., v1.0)"
  echo "-y path_to_rufus_container   If not provided, will look in current directory for rufus.sif"
  echo "-z rufus_threads	Number of threads provided to RUFUS; defaults to 36 for entire genome; 10 for 1MB windows (NOTE: must be less than cpus_per_call)"
  echo "-e email  The email address to notify with slurm updates"
  echo "-q slurm_job_queue_limit    The maximum amount of jobs able to be ran at once; defaults to 20"
  echo "-t slurm_time_limit   The maximum amount of time to let the slurm job run; defaults to 7 days for full run, or one hour per window (DD-HH:MM:SS)"
  echo "-M memory_per_call    How much memory to allot to the rufus calling stage job; default 150G for entire genome; 20G for 1MB windows (e.g. 150G or 20G)"
  echo "-C cpus_per_call      How many cpus to allot to each rufus calling stage job; default 40 for entire genome; 12 for 1MB windows"
  echo "-h help	Print usage"
  echo ""
  echo "Output files are written to the current working directory."
	exit 1
}

# Initialize variables
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
REFERENCE_HASH_RUFUS_ARG=""
KG1_HASH_DIR=""
KG1_HASH_VERSION=""
CONTROL_HASH_DIR=""
CONTROL_HASH_VERSION=""
MEM_PER_JOB=""
CPUS_PER_JOB=""

# Parse command line options using getopts
while getopts ":s:c:b:a:p:r:m:w:e:l:q:t:f:x:y:z:h:M:CK:G:D:V:" opt; do
    case ${opt} in
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
        K)
            KG1_HASH_DIR=$OPTARG
            ;;
        G)
            KG1_HASH_VERSION=$OPTARG
            ;;
        D)
            CONTROL_HASH_DIR=$OPTARG
            ;;
        V)
            CONTROL_HASH_VERSION=$OPTARG
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
if [[ -z "$SUBJECT_RUFUS_ARG" || -z "$GENOME_BUILD_RUFUS_ARG" || -z "$REFERENCE_RUFUS_ARG" || -z "$SLURM_ACCOUNT_RUFUS_ARG" || -z "$SLURM_PARTITION_RUFUS_ARG" || -z "$SLURM_ARRAY_JOB_LIMIT_RUFUS_ARG" ]]; then
    echo "ERROR: Missing required argument(s); please see usage instructions with -h." >&2
    exit 1
fi

# Check that subject file exists
if [ ! -f "$SUBJECT_RUFUS_ARG" ]; then
	echo "ERROR: subject file $SUBJECT_RUFUS_ARG does not exist or cannot be read." >&2
	exit 1
fi

# Check that all control files exist
for control in "${CONTROLS_RUFUS_ARG[@]}"; do
	if [ ! -f "$control" ]; then
		echo "ERROR: control file $control does not exist or cannot be read." >&2
		exit 1
	else
		if [ "$CONTROL_STRING_RUFUS_ARG" == "" ]; then
			CONTROL_STRING_RUFUS_ARG="$control"
		else
			CONTROL_STRING_RUFUS_ARG="${CONTROL_STRING_RUFUS_ARG}, $control"
		fi
	fi
done

# Validate that subject and control files are all the same type (bam, cram, or fastq)
get_input_type() {
	case "$1" in
		*.cram) echo "cram" ;;
		*.bam) echo "bam" ;;
		*.fastq.gz|*.fq.gz) echo "fastq" ;;
		*.fastq|*.fq) echo "fastq" ;;
		*) echo "unknown" ;;
	esac
}

SUBJECT_TYPE=$(get_input_type "$SUBJECT_RUFUS_ARG")
if [ "$SUBJECT_TYPE" == "unknown" ]; then
	echo "ERROR: subject file $SUBJECT_RUFUS_ARG has an unrecognized file type. Supported types: .bam, .cram, .fastq, .fq, .fastq.gz, .fq.gz" >&2
	exit 1
fi
for control in "${CONTROLS_RUFUS_ARG[@]}"; do
	CTRL_TYPE=$(get_input_type "$control")
	if [ "$CTRL_TYPE" != "$SUBJECT_TYPE" ]; then
		echo "ERROR: all subject and control files must be the same type, but subject is ${SUBJECT_TYPE} and control $control is ${CTRL_TYPE}. Please ensure all inputs are either all BAMs, all CRAMs, or all FASTQs." >&2
		exit 1
	fi
done

# Check that reference file exists
if [ ! -f "$REFERENCE_RUFUS_ARG" ]; then
	echo "ERROR: reference file $REFERENCE_RUFUS_ARG does not exist or cannot be read." >&2
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

if [ "$THREAD_LIMIT_RUFUS_ARG" -ge "$CPUS_PER_JOB" ]; then
	echo "ERROR: thread limit ($THREAD_LIMIT_RUFUS_ARG) must be less than cpus per job ($CPUS_PER_JOB)." >&2
	exit 1
fi
# Validate per-region hash flags: cannot specify both local dir and S3 version for same type
if [ -n "$KG1_HASH_DIR" ] && [ -n "$KG1_HASH_VERSION" ]; then
    echo "ERROR: Cannot specify both -K (local KG1 hash dir) and -G (KG1 S3 version). Use one or the other." >&2
    exit 1
fi
if [ -n "$CONTROL_HASH_DIR" ] && [ -n "$CONTROL_HASH_VERSION" ]; then
    echo "ERROR: Cannot specify both -D (local control hash dir) and -V (control S3 version). Use one or the other." >&2
    exit 1
fi

# Validate local hash directories exist if provided
if [ -n "$KG1_HASH_DIR" ] && [ ! -d "$KG1_HASH_DIR" ]; then
    echo "ERROR: KG1 hash directory does not exist: $KG1_HASH_DIR" >&2
    exit 1
fi
if [ -n "$CONTROL_HASH_DIR" ] && [ ! -d "$CONTROL_HASH_DIR" ]; then
    echo "ERROR: Control hash directory does not exist: $CONTROL_HASH_DIR" >&2
    exit 1
fi

# Check that if path to image not provided, it's in the current dir
if [ -z $CONTAINER_PATH_RUFUS_ARG ]; then
	if [ ! -f "rufus.sif" ]; then
		echo "ERROR: rufus.sif not in current directory - please provide path to container or put it in this one under rufus.sif"
	fi
fi

# Collect unique parent directories for all input files to use as bind mounts.
# Singularity --bind preserves host paths inside the container (no remapping needed).
collect_bind_dirs() {
    local -A seen_dirs
    local dirs=()

    # Always include pwd for output
    seen_dirs["$(pwd)"]=1
    dirs+=("$(pwd)")

    local files=("$SUBJECT_RUFUS_ARG" "$REFERENCE_RUFUS_ARG")
    for control in "${CONTROLS_RUFUS_ARG[@]}"; do
        files+=("$control")
    done
    if [ -n "$REFERENCE_HASH_RUFUS_ARG" ]; then
        files+=("$REFERENCE_HASH_RUFUS_ARG")
    fi
    for exclude in "${EXCLUDE_HASH_LIST_RUFUS_ARG[@]}"; do
        files+=("$exclude")
    done

    # Add hash directories directly (not individual files)
    local hash_dirs=()
    if [ -n "$KG1_HASH_DIR" ]; then
        hash_dirs+=("$KG1_HASH_DIR")
    fi
    if [ -n "$CONTROL_HASH_DIR" ]; then
        hash_dirs+=("$CONTROL_HASH_DIR")
    fi

    for hd in "${hash_dirs[@]}"; do
        local resolved_hd
        resolved_hd="$(realpath "$hd")"
        if [ -z "${seen_dirs[$resolved_hd]+x}" ]; then
            seen_dirs["$resolved_hd"]=1
            dirs+=("$resolved_hd")
        fi
    done

    for f in "${files[@]}"; do
        local d
        d="$(dirname "$(realpath "$f")")"
        if [ -z "${seen_dirs[$d]+x}" ]; then
            seen_dirs["$d"]=1
            dirs+=("$d")
        fi
    done

    # Join with commas
    local IFS=','
    echo "${dirs[*]}"
}

BIND_MOUNTS="$(collect_bind_dirs)"

# Export variables for use in the main script
export BIND_MOUNTS
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
export REFERENCE_HASH_RUFUS_ARG
export KG1_HASH_DIR
export KG1_HASH_VERSION
export CONTROL_HASH_DIR
export CONTROL_HASH_VERSION
export MEM_PER_JOB
export CPUS_PER_JOB