#!/bin/bash

# CONSTANTS
PR_WORKER="/opt/RUFUS/aws_launch/process_region_worker.sh"
DEV_MOUNT="-v /home/ubuntu/RUFUS/runRufus.sh:/opt/RUFUS/runRufus.sh \
  -v /home/ubuntu/RUFUS/scripts:/opt/RUFUS/scripts \
  -v /home/ubuntu/RUFUS/resource_helpers:/opt/RUFUS/resource_helpers \
  -v /home/ubuntu/RUFUS/post_process:/opt/RUFUS/post_process \
  -v /home/ubuntu/RUFUS/resources:/opt/RUFUS/resources
  -v /home/ubuntu/RUFUS/bin/RUFUS.interpret:/opt/RUFUS/bin/RUFUS.interpret"
#DEV_MOUNT=""

# Check for required argument
ENV_FILE="$1"
if [ -z "$ENV_FILE" ]; then
    echo "Error: Please provide PATH_TO_RUFUS_ENV argument" >&2
    exit 1
fi

# Check RUFUS env file arg actually exists
if [ -f "$ENV_FILE" ]; then
    set -a
    source <(grep -v '^#' $ENV_FILE | grep -v '^[[:space:]]*$' | sed 's/\r$//')
    set +a
else
    echo "Error: $ENV_FILE file not found - please provide valid path to rufus.env file"
    exit 1
fi

# Checks for required arguments and formatting + bounds of integer arguments
check_inputs() {
    # Check for all required variables to be filled in rufus.env
    REQUIRED_VARS=(SUBJECT_FILE KMER_DEPTH_CUTOFF THREAD_LIMIT REFERENCE_FASTA RUFUS_DOCKER_IMAGE)
    for var in "${REQUIRED_VARS[@]}"; do
        if [ -z "${!var}" ]; then
            echo "Error: Required variable $var is not set in rufus.env"
            exit 1
        fi
    done

    # Check for subject file existence
    subject_path="${SUBJECT_FILE}"
    if [ ! -f "$subject_path" ]; then
        echo "Error: SUBJECT_FILE $subject_path not found" >&2
        exit 1
    fi

    # Check for control files existence
    CONTROL_FILE_ARRAY=()
    for control in "${CONTROLS[@]}"; do
        if [ ! -f "$control" ]; then
            echo "Error: Control file $control not found" >&2
            exit 1
        else
            CONTROL_FILE_ARRAY+=("$control")
        fi
    done

    # Check for reference fasta existence
    reference_path="${REFERENCE_FASTA}"
    if [ ! -f "$reference_path" ]; then
        echo "Error: REFERENCE_FASTA $reference_path not found" >&2
        exit 1
    fi

    # Check if region file provided, that it exists
    if [ -n "$REGION_FILE" ]; then
        if [ ! -f "$REGION_FILE" ]; then
            echo "Error: REGION_FILE $REGION_FILE not found. Please provide valid file or leave empty for whole genome mode." >&2
            exit 1
        else
            REGION_PATH="$REGION_FILE"
        fi
    fi

    # Check that thread limit is an int greater than 0
    if ! [[ "$THREAD_LIMIT" =~ ^[0-9]+$ ]] || [ "$THREAD_LIMIT" -le 0 ]; then
        echo "Error: THREAD_LIMIT must be a positive integer" >&2
        exit 1
    fi

    # Check that if job threshold is present, it is an int greater than 0
    if [ -n "$JOB_THRESHOLD" ]; then
        if ! [[ "$JOB_THRESHOLD" =~ ^[0-9]+$ ]] || [ "$JOB_THRESHOLD" -le 0 ]; then
            echo "Error: JOB_THRESHOLD must be a positive integer" >&2
            exit 1
        fi
    else
        JOB_THRESHOLD=1
    fi

    # Check that kmer depth cutoff is an int and warn if less than 3
    if ! [[ "$KMER_DEPTH_CUTOFF" =~ ^[0-9]+$ ]]; then
        echo "Error: KMER_DEPTH_CUTOFF must be a positive integer" >&2
        exit 1
    elif [ "$KMER_DEPTH_CUTOFF" -lt 3 ]; then
        echo "Warning: KMER_DEPTH_CUTOFF is set to less than 3, which may lead to increased false positives" >&2
    fi

    # Check if KMER_LENGTH filled out, is an int and warn if not 25
    if [ -n "$KMER_LENGTH" ]; then
        if ! [[ "$KMER_LENGTH" =~ ^[0-9]+$ ]]; then
            echo "Error: KMER_LENGTH must be a positive integer" >&2
            exit 1
        elif [ "$KMER_LENGTH" -ne 25 ]; then
            echo "Warning: KMER_LENGTH is set to $KMER_LENGTH, RUFUS has been robustly tested with a KMER_LENGTH of 25 and is recommended" >&2
        fi
    fi

    # Check if WINDOW_SIZE filled out, if it is, make sure an int and is 1000
    if [ -n "$WINDOW_SIZE" ]; then
        if ! [[ "$WINDOW_SIZE" =~ ^[0-9]+$ ]]; then
            echo "Error: WINDOW_SIZE must be a positive integer" >&2
            exit 1
        elif [ "$WINDOW_SIZE" -ne 1000 ]; then
            echo "ERROR: WINDOW_SIZE is set to $WINDOW_SIZE - must be either 1000 or empty" >&2
            exit 1
        fi
    fi

    # If we don't have a working dir, set it to .
    if [ -z "$WORKING_DIR" ]; then
        WORKING_DIR="."
    fi
}
export -f check_inputs

# Check for correct controls setup and returns paths needed for mounting if necessary
set_up_controls() {
    # Check for controls here and notify if using internal
    if [ ${#CONTROLS[@]} -eq 0 ]; then
        echo "No paired controls provided, running RUFUS in internal control mode..." >&2
    fi

    local mount_clause=""
    # Make sure the local directory exists if provided
    if [ ! -z "${CONTROL_HASH_LOCAL_DIR}" ]; then
        if [ ! -d "${CONTROL_HASH_LOCAL_DIR}" ]; then
            echo "Error: CONTROL_HASH_LOCAL_DIR ${CONTROL_HASH_LOCAL_DIR} does not exist. Please ensure directory exists or leave CONTROL_HASH_LOCAL_DIR empty for S3 fetching." >&2
            return 1
        else
            # Make sure directory has at least one *.Jhash file in it
            shopt -s nullglob
            jhash_files=("${CONTROL_HASH_LOCAL_DIR}"/*.Jhash)
            shopt -u nullglob
            if [ ${#jhash_files[@]} -eq 0 ]; then
                echo "Error: CONTROL_HASH_LOCAL_DIR ${CONTROL_HASH_LOCAL_DIR} does not contain any *.Jhash files. Please ensure directory has Jhash files or leave CONTROL_HASH_LOCAL_DIR empty for S3 fetching." >&2
                return 1
            fi
            mount_clause="-v ${CONTROL_HASH_LOCAL_DIR}:/mnt/rufus_resources/control_hashes "
        fi
    fi

    # If we have paired controls provided, also use those
    if [ "${#CONTROL_FILE_ARRAY[@]}" -ne 0 ]; then
        # Concatenate controls into -c delimited string
        ctrl_arg=""
        for control in "${CONTROL_FILE_ARRAY[@]}"; do
            ctrl_arg+="-v /mnt/$control "
        done
        mount_clause+="$ctrl_arg"
    fi
    echo "$mount_clause"
}
export -f set_up_controls

# Check for correct 1000G setup and returns paths needed for mounting if necessary
set_up_kg1() {
    local mount_clause=""

    # Notify if not removing 1000G variants
    if [ "$NO_KG1_REMOVAL" == "true" ] || [ "$NO_KG1_REMOVAL" == "TRUE" ]; then
        echo "Warning: not removing common population variants in the 1000G cohort" >&2
    else
        # Make sure the local directory exists if provided
        if [ ! -z "${KG1_HASH_LOCAL_DIR}" ]; then
            if [ ! -d "${KG1_HASH_LOCAL_DIR}" ]; then
                echo "Error: KG1_HASH_LOCAL_DIR ${KG1_HASH_LOCAL_DIR} does not exist. Please ensure directory exists or leave KG1_HASH_LOCAL_DIR empty for S3 fetching." >&2
                return 1
            else
                # Make sure directory has at least one *.Jhash file in it
                shopt -s nullglob
                jhash_files=("${KG1_HASH_LOCAL_DIR}"/*.Jhash)
                shopt -u nullglob
                if [ ${#jhash_files[@]} -eq 0 ]; then
                    echo "Error: KG1_HASH_LOCAL_DIR ${KG1_HASH_LOCAL_DIR} does not contain any *.Jhash files. Please ensure directory has Jhash files or leave KG1_HASH_LOCAL_DIR empty for S3 fetching." >&2
                    return 1
                fi
                mount_clause="-v ${KG1_HASH_LOCAL_DIR}:/mnt/rufus_resources/kg1_hashes"
            fi
        fi
    fi
    echo "$mount_clause"
}
export -f set_up_kg1

set_up_ref() {
    local ref_base=$(basename ${REFERENCE_FASTA})
    local mount_clause="-v ${REFERENCE_FASTA}:/mnt/${ref_base} "
    local build_refs="FALSE"

    # Determine the base name to check once
    if [[ "$REFERENCE_FASTA" == *.gz ]]; then
        ref_file="${REFERENCE_FASTA%.gz}"
    else
        ref_file="$REFERENCE_FASTA"
    fi

    # Check all required index files in one loop
    index_mounts=""
    for ext in sa bwt pac amb ann; do
        if [[ ! -e "${ref_file}.${ext}" ]]; then
            build_refs="TRUE"
            break
        else
            ref_base=$(basename ${ref_file})
            index_mounts+="-v ${ref_file}.${ext}:/mnt/${ref_base}.${ext} "
        fi
    done

    # Mount indexes only if we have them all and don't need to build
    if [ "$build_refs" == "FALSE" ]; then
        mount_clause+="$index_mounts"
    fi

    echo "$mount_clause|$build_refs"
}
export -f set_up_ref

# Start work
check_inputs
IFS='|' read -r ref_mount build_refs < <(set_up_ref)
control_mount=$(set_up_controls) || exit 1
kg1_mount=$(set_up_kg1) || exit 1
subject_base=$(basename ${SUBJECT_FILE})
input_mount_clause="$ref_mount $control_mount $kg1_mount -v ${SUBJECT_FILE}:/mnt/${subject_base}"

CONTAINER_ID=$(docker run -d --rm --name rufus-worker \
  -v ${WORKING_DIR}:/mnt \
  -v ${ENV_FILE}:/mnt/rufus_resources/rufus.env \
  $input_mount_clause \
  $DEV_MOUNT \
  $RUFUS_DOCKER_IMAGE \
  tail -f /dev/null)

start_time=$(date +%s)

# Write commands for final vcf (do inside container so have access to RUFUS versioning)
docker exec ${CONTAINER_ID} bash /opt/RUFUS/resource_helpers/write_command_args.sh "$CONTAINER_ID"

# Check for BWA indexes and create if necessary
if [ "$build_refs" == "TRUE" ]; then
    echo "Generating BWA indexes for reference fasta..."
    docker exec ${CONTAINER_ID} bash /opt/RUFUS/resource_helpers/check_for_bwa_indexes.sh "${REFERENCE_FASTA}"
fi

# Make rufus resource dirs inside container
docker exec ${CONTAINER_ID} mkdir -p /mnt/rufus_resources/control_hashes
docker exec ${CONTAINER_ID} mkdir -p /mnt/rufus_resources/kg1_hashes
docker exec ${CONTAINER_ID} mkdir -p /mnt/rufus_resources/logs

# Start work
echo "Starting RUFUS job(s)..."
if [ -n "$REGION_FILE" ]; then
    parallel -j "$JOB_THRESHOLD" "${PR_WORKER}" "$CONTAINER_ID" {} :::: "$REGION_FILE"
else
    "${PR_WORKER}" "$CONTAINER_ID" ""
fi

concat_ctrl_post_arg=""
if [ ${#CONTROLS[@]} -gt 0 ]; then
    # Concatenate controls without -c delimiters
    concat_ctrls=""
    for control in "${CONTROL_FILE_ARRAY[@]}"; do
        ctrl_base=$(basename ${control})
        concat_ctrls+="/mnt/$ctrl_base "
    done
    concat_ctrl_post_arg="-c $concat_ctrls"
fi

subject_base=$(basename ${SUBJECT_FILE})
ref_base=$(basename ${REFERENCE_FASTA})

# Wait for all jobs to finish before combining + post-processing
echo "All RUFUS regional jobs completed. Starting merge and post-process..."
docker exec ${CONTAINER_ID} bash /opt/RUFUS/post_process/post_process.sh -s "/mnt/$subject_base" -r "/mnt/${ref_base}" -w "$WINDOW_SIZE" -d "/mnt" "$concat_ctrl_post_arg"

# Stop container and clean up
docker stop rufus-worker

end_time=$(date +%s)
elapsed=$((end_time - start_time))
echo "RUFUS completed. Total run time: $elapsed"