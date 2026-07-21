#!/bin/bash

# Statics
DEFAULT_1MB_CPUS_PER_JOB="12"
DEFAULT_1MB_MEM_PER_JOB="20G"
DEFAULT_WG_CPUS_PER_JOB="40"
DEFAULT_WG_MEM_PER_JOB="150G"


usage() {
  echo "Usage: $0 [-s subject1,subject2,...] [-c control1,control2,control3...] [-b genome_build] [-a slurm_account] [-p slurm_partition] ...options"
  echo "Required Arguments:"
  echo "-s subject(s) A single subject or comma-delimited array of multiple subject BAM/CRAM files (full paths)"
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
  echo "-d dev_binds     Comma-delimited list of host:container bind mounts for dev testing (e.g., /local/runRufus.sh:/opt/RUFUS/runRufus.sh)"
  echo "-P par_low_cov_threshold  Control k-mer count ceiling below which a variant is flagged as low-coverage-parent/inherited (default 7; set to 0 to disable, e.g. when using an assembly as the control)"
  echo "-h help	Print usage"
  echo ""
  echo "Output files are written to the current working directory."
	exit 1
}

# Initialize variables
SUBJECTS_RUFUS_ARG=()
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
PAR_LOW_COV_THRESHOLD_RUFUS_ARG="7"
DEV_BIND_MOUNTS_ARG=()

# Parse command line options using getopts
while getopts ":s:c:b:a:p:r:m:w:e:l:q:t:f:x:y:z:h:M:CK:G:D:V:d:P:" opt; do
    case ${opt} in
        s)
            IFS=',' read -r -a SUBJECTS_RUFUS_ARG <<< "$OPTARG"
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
        d)
            IFS=',' read -r -a DEV_BIND_MOUNTS_ARG <<< "$OPTARG"
            ;;
        P)
            if ! [[ "$OPTARG" =~ ^[0-9]+$ ]]; then
                echo "ERROR: -P par_low_cov_threshold must be a non-negative integer." >&2
                exit 1
            fi
            PAR_LOW_COV_THRESHOLD_RUFUS_ARG=$OPTARG
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
if [[ ${#SUBJECTS_RUFUS_ARG[@]} -eq 0 || -z "$GENOME_BUILD_RUFUS_ARG" || -z "$REFERENCE_RUFUS_ARG" || -z "$SLURM_ACCOUNT_RUFUS_ARG" || -z "$SLURM_PARTITION_RUFUS_ARG" || -z "$SLURM_ARRAY_JOB_LIMIT_RUFUS_ARG" ]]; then
    echo "ERROR: Missing required argument(s); please see usage instructions with -h." >&2
    exit 1
fi

# Check that all subject files exist
for subject in "${SUBJECTS_RUFUS_ARG[@]}"; do
	if [ ! -f "$subject" ]; then
		echo "ERROR: subject file $subject does not exist or cannot be read." >&2
		exit 1
	fi
done

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

# Validate type of first subject, then ensure all subjects + controls match
SUBJECT_TYPE=$(get_input_type "${SUBJECTS_RUFUS_ARG[0]}")
if [ "$SUBJECT_TYPE" == "unknown" ]; then
	echo "ERROR: subject file ${SUBJECTS_RUFUS_ARG[0]} has an unrecognized file type. Supported types: .bam, .cram, .fastq, .fq, .fastq.gz, .fq.gz" >&2
	exit 1
fi
for subject in "${SUBJECTS_RUFUS_ARG[@]:1}"; do
	SUBJ_TYPE=$(get_input_type "$subject")
	if [ "$SUBJ_TYPE" != "$SUBJECT_TYPE" ]; then
		echo "ERROR: all subject files must be the same type, but first subject is ${SUBJECT_TYPE} and $subject is ${SUBJ_TYPE}. Please ensure all inputs are either all BAMs, all CRAMs, or all FASTQs." >&2
		exit 1
	fi
done
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

# Check that BWA indexes exist alongside the reference.
# runRufus.sh prefers the extension-stripped prefix when <prefix>.sa is present and
# otherwise falls back to the reference path itself, so validate whichever it will pick.
# Without this, the missing index only surfaces at the bwa mem step, which in region
# mode is hours into every queued array task.
REFERENCE_BWA_PREFIX="$REFERENCE_RUFUS_ARG"
if [ -e "${REFERENCE_RUFUS_ARG%.*}.sa" ]; then
	REFERENCE_BWA_PREFIX="${REFERENCE_RUFUS_ARG%.*}"
fi

MISSING_BWA_INDEXES=()
for bwa_suffix in amb ann bwt pac sa; do
	if [ ! -e "${REFERENCE_BWA_PREFIX}.${bwa_suffix}" ]; then
		MISSING_BWA_INDEXES+=("${REFERENCE_BWA_PREFIX}.${bwa_suffix}")
	fi
done
if [ ! -e "${REFERENCE_RUFUS_ARG}.fai" ]; then
	MISSING_BWA_INDEXES+=("${REFERENCE_RUFUS_ARG}.fai")
fi

if [ ${#MISSING_BWA_INDEXES[@]} -ne 0 ]; then
	echo "ERROR: reference $REFERENCE_RUFUS_ARG is missing required index files:" >&2
	for missing in "${MISSING_BWA_INDEXES[@]}"; do
		echo "       $missing" >&2
	done
	echo "RUFUS aligns candidate reads with BWA and cannot run without these." >&2
	echo "Build them once with:" >&2
	echo "       bash \${RUFUS_ROOT}/resource_helpers/build_bwa_indexes.sh $REFERENCE_RUFUS_ARG" >&2
	echo "(indexing a human-sized reference takes roughly an hour)" >&2
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
	if [ -z "$SLURM_TIME_LIMIT_RUFUS_ARG" ]; then
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

# Validate dev bind mounts and build singularity --bind args string
DEV_BIND_ARGS=""
if [ ${#DEV_BIND_MOUNTS_ARG[@]} -gt 0 ]; then
    for bind_spec in "${DEV_BIND_MOUNTS_ARG[@]}"; do
        host_path="${bind_spec%%:*}"
        container_path="${bind_spec#*:}"
        if [ "$host_path" == "$bind_spec" ]; then
            echo "ERROR: dev bind mount '$bind_spec' must be in host:container format (e.g., /local/runRufus.sh:/opt/RUFUS/runRufus.sh)" >&2
            exit 1
        fi
        if [ ! -e "$host_path" ]; then
            echo "ERROR: dev bind mount host path does not exist: $host_path" >&2
            exit 1
        fi
        DEV_BIND_ARGS+=" --bind ${bind_spec}"
    done
    echo "DEV MODE: additional bind mounts:${DEV_BIND_ARGS}"
fi

# Collect unique parent directories for all input files to use as bind mounts.
# Singularity --bind preserves host paths inside the container (no remapping needed).
collect_bind_dirs() {
    local -A seen_dirs
    local dirs=()

    # Always include pwd for output
    seen_dirs["$(pwd)"]=1
    dirs+=("$(pwd)")

    local files=("${SUBJECTS_RUFUS_ARG[@]}" "$REFERENCE_RUFUS_ARG")
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
export SUBJECTS_RUFUS_ARG
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
export PAR_LOW_COV_THRESHOLD_RUFUS_ARG
export DEV_BIND_ARGS