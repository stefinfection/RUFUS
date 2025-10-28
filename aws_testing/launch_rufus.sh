#!/bin/bash

# Check for RUFUS env file
if [ -f aws_rufus.env ]; then
    set -a
    source <(grep -v '^#' aws_rufus.env | grep -v '^[[:space:]]*$' | sed 's/\r$//')
    set +a
else
    echo "Error: rufus.env file not found"
    exit 1
fi

# Check for all required variables
REQUIRED_VARS=(HOST_DATA_DIR CONTAINER_PATH SUBJECT_FILE REGION_FILE KMER_LENGTH KMER_DEPTH_CUTOFF THREAD_LIMIT JOB_THRESHOLD OTHER_FLAGS REFERENCE_FASTA CONTROL_HASH_VERSION KG1_HASH_VERSION WINDOW_SIZE)
for var in "${REQUIRED_VARS[@]}"; do
    if [ -z "${!var}" ]; then
        echo "Error: Required variable $var is not set in rufus.env"
        exit 1
    fi
done

# Helper functions
get_control_hash() {
    local region=$1
    
    # Check for resources directory
    if [ ! -d "${HOST_DATA_DIR}/rufus_resources/control_hashes" ]; then
        echo -n "Error: Internal control resources directory not found at ${HOST_DATA_DIR}/rufus_resources/control_hashes - "
        echo "Please ensure the RUFUS resources directory is copied or soft-linked within the HOST_DATA_DIR assigned in the rufus.env file"
        exit 1
    fi
    
    local ctrl_hash=""
    if [ "$region" == "" ]; then
        # If we don't have a region, use entire genome wide Jhash
        ctrl_hash="/mnt/rufus_resources/control_hashes/wg.control_${CONTROL_HASH_VERSION}.Jhash"
    else
        # Otherwise, use region specific technical control
        hash_arg=$(echo "$region" | tr ':-' '_')
        ctrl_hash="/mnt/rufus_resources/control_hashes/${hash_arg}.control_${CONTROL_HASH_VERSION}.Jhash"
    fi

    echo "$ctrl_hash"
}
export -f get_internal_control

get_kg1_hash() {
    local region=$1

    local kg1_hash_arg=""
    # Check for resources directory
    if [ ! -d "${HOST_DATA_DIR}/rufus_resources/kg1_hashes" ]; then
        echo -n "Error - rufus_resources/kg1_hashes directory not found in ${HOST_DATA_DIR} - "
        echo "Please ensure the RUFUS resources directory is copied or soft-linked within the HOST_DATA_DIR assigned in the rufus.env file"
        exit 1
    fi

    local hash_arg=""
    if [ "$region" == "" ]; then
        # If we don't have a region, use entire genome wide Jhash
        kg1_hash_arg="/mnt/rufus_resources/kg1_hashes/wg.${KG1_HASH_VERSION}.Jhash"
    else
        # Convert chrN:n-m to chrN_n_m
        kg1_hash_arg=$(echo "$region" | tr ':-' '_')
        kg1_hash_arg="/mnt/rufus_resources/kg1_hashes/${hash_arg}.${KG1_HASH_VERSION}.Jhash"
    fi

    echo "$kg1_hash_arg"
}
export -f get_kg1_hash

get_reference() {
    if [ ! -d "${HOST_DATA_DIR}/rufus_resources/references" ] || [ ! -f "${HOST_DATA_DIR}/rufus_resources/references/${REFERENCE_FASTA}" ]; then
        # Check for resources directory and that the reference is in there - warn if not
        echo -n "Warning - ${REFERENCE_FASTA} not found in ${HOST_DATA_DIR}/rufus_resources/references directory - "
        echo "Utilizing pre-built BWA references provided in the RUFUS resources directory can speed up run time next time"
        reference="mnt/${REFERENCE_FASTA}"
    else
        # If not, just use what they provided
        reference="/mnt/rufus_resources/references/${REFERENCE_FASTA}"
    fi

    echo "$reference"
}

process_region() {
    local region=$1

    ctrl_arg=""
    # Check to see if controls are provided
    if [ "${#CONTROL_FILE_ARRAY[@]}" -eq 0 ]; then
        echo "No control samples provided, using internal control for single sample mode"
        internal_ctrl=$($get_control_hash $region)
        ctrl_arg="-c $internal_ctrl "
    else
        # Concatenate controls into a single -c delimited string
        for control in "${CONTROL_FILE_ARRAY[@]}"; do
            ctrl_arg+="-c /mnt/$control "
        done
    fi

    kg1_hash=$(get_kg1_hash $region)
    kg1_hash_arg="-e $kg1_hash"

    ref=$(get_reference)
    ref_arg="-r $ref"

    echo "singularity exec --bind ${HOST_DATA_DIR}:/mnt ${CONTAINER_PATH} bash /opt/RUFUS/runRufus.sh \
        -s /mnt/$SUBJECT_FILE \
        $ctrl_arg \
        $ref_arg \
        -m $KMER_DEPTH_CUTOFF \
        -k $KMER_LENGTH \
        -t $THREAD_LIMIT \
        $OTHER_FLAGS \
        $kg1_hash_arg \
        $region"
}
export -f process_region

# Check for regions file to determine if running in WG mode or region mode
if [ -e "$REGION_FILE" ]; then
    num=$(cat "$REGION_FILE" | wc -l)
    first=$(cat "$REGION_FILE" | head -n 1)
    last=$(cat "$REGION_FILE" | tail -n 1)
    echo "Using $REGION_FILE with $num regions - first is $first and last is $last"
else
    touch empty_region.txt
    REGION_FILE="empty_region.txt"
    echo "Running RUFUS in whole genome mode"
fi

# TODO: port to Docker command
# Start work
echo "Starting parallel RUFUS jobs..."
parallel -j "$JOB_THRESHOLD" process_region {} :::: "$REGION_FILE"

# Concatenate controls without -c delimiters
concat_ctrl_post_arg=""
for control in "${CONTROL_FILE_ARRAY[@]}"; do
    concat_ctrl_post_arg+="/mnt/$control "
done

ref=$($get_reference)

# Wait for all jobs to finish before combining + post-processing
#echo "All RUFUS regional jobs completed. Starting merge and post-process..."
#singularity exec --bind "${HOST_DATA_DIR}:/mnt" "${CONTAINER_PATH}$" bash /opt/RUFUS/post_process/post_process.sh -s "/mnt/$SUBJECT_FILE" -c "$concat_ctrl_post_arg" -r "$ref" -w "$WINDOW_SIZE" -d "/mnt"
echo "RUFUS completed. Total run time: "