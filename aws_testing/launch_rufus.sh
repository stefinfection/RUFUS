#!/bin/bash

$ENV_FILE=$1

# Check required arg
if [ -z "$ENV_FILE" ]; then
    echo "Usage: $0 /path/to/rufus.env"
    exit 1
fi

start_time=$(date +%s)

# Check for RUFUS env file arg actually exists
if [ -f "$ENV_FILE" ]; then
    set -a
    source <(grep -v '^#' $ENV_FILE | grep -v '^[[:space:]]*$' | sed 's/\r$//')
    set +a
else
    echo "Error: $ENV_FILE file not found"
    exit 1
fi

# Check for all required variables to be filled in rufus.env
REQUIRED_VARS=(HOST_DATA_DIR CONTAINER_PATH SUBJECT_FILE REGION_FILE KMER_LENGTH KMER_DEPTH_CUTOFF THREAD_LIMIT JOB_THRESHOLD OTHER_FLAGS REFERENCE_FASTA CONTROL_HASH_VERSION KG1_HASH_VERSION WINDOW_SIZE)
for var in "${REQUIRED_VARS[@]}"; do
    if [ -z "${!var}" ]; then
        echo "Error: Required variable $var is not set in rufus.env"
        exit 1
    fi
done

# Fetch functions to get resources from S3
fetch_kg1_hash() {
    local region=$1
    fmtd_reg=$(echo "$region" | tr ':-' '_')

    if [ ! -d "${HOST_DATA_DIR}/rufus_resources" ]; then
        mkdir -p "${HOST_DATA_DIR}/rufus_resources"
    fi

    local kg1_hash=""
    if [ "$region" == "" ]; then
        # If we don't have a region, use entire genome wide Jhash
        aws s3 sync "s3://rufus.marth.lab/public_access_data/kg1_hashes/${KG1_HASH_VERSION}/wg_kg1_${KG1_HASH_VERSION}.Jhash" "${HOST_DATA_DIR}/rufus_resources/wg_kg1_${KG1_HASH_VERSION}.Jhash"
        kg1_hash="mnt/rufus_resources/wg_kg1_${KG1_HASH_VERSION}.Jhash"
    else
        # Convert chrN:n-m to chrN_n_m
        chrom=$(echo "$fmtd_reg" | cut -d'_' -f1)
        aws s3 sync "s3://rufus.marth.lab/public_access_data/kg1_hashes/${KG1_HASH_VERSION}/${chrom}/${fmtd_reg}_kg1_${KG1_HASH_VERSION}.Jhash" "${HOST_DATA_DIR}/rufus_resources/${fmtd_reg}_kg1_${KG1_HASH_VERSION}.Jhash"
        kg1_hash="/mnt/rufus_resources/${fmtd_reg}_kg1_${KG1_HASH_VERSION}.Jhash"
    fi

    echo "$kg1_hash"
}

fetch_control_hash() {
    local region=$1
    fmtd_reg=$(echo "$region" | tr ':-' '_')

    if [ ! -d "${HOST_DATA_DIR}/rufus_resources" ]; then
        mkdir -p "${HOST_DATA_DIR}/rufus_resources"
    fi

    local ctrl_hash=""
    if [ "$region" == "" ]; then
        # If we don't have a region, use entire genome wide Jhash
        aws s3 sync "s3://rufus.marth.lab/public_access_data/control_hashes/${CONTROL_HASH_VERSION}/wg_control_${CONTROL_HASH_VERSION}.Jhash" "${HOST_DATA_DIR}/rufus_resources/wg_control_${CONTROL_HASH_VERSION}.Jhash"
        ctrl_hash="mnt/rufus_resources/wg_control_${CONTROL_HASH_VERSION}.Jhash"
    else
        # Convert chrN:n-m to chrN_n_m
        chrom=$(echo "$fmtd_reg" | cut -d'_' -f1)
        aws s3 sync "s3://rufus.marth.lab/public_access_data/kg1_hashes/${CONTROL_HASH_VERSION}/${chrom}/${fmtd_reg}_control_${CONTROL_HASH_VERSION}.Jhash" "${HOST_DATA_DIR}/rufus_resources/${fmtd_reg}_control_${CONTROL_HASH_VERSION}.Jhash"
        ctrl_hash="/mnt/rufus_resources/${fmtd_reg}_control_${CONTROL_HASH_VERSION}.Jhash"
    fi
    echo "$ctrl_hash"
}

# Helper functions to get resources locally

# Checks to see if we already have BWA indexes premade for the reference argument, and points there if so
get_reference() {
    reference="${HOST_DATA_DIR}/${REFERENCE_FASTA}" 

    # If we have the exact reference BWA indexes already, point there to save some time
    if [ -d "${HOST_DATA_DIR}/rufus_resources/references" ] && [ -f "${HOST_DATA_DIR}/rufus_resources/references/${REFERENCE_FASTA}" ]; then
        reference="/mnt/rufus_resources/references/${REFERENCE_FASTA}"
    fi

    echo "$reference"
}

# Looks for 1000G hashes in ${HOST_DATA_DIR}/rufus_resources/kg1_hashes
# If can't find, will pull from S3
get_kg1_hash() {
    local region=$1
    local kg1_hash_arg=""

    if [ ! -d "${HOST_DATA_DIR}/rufus_resources/kg1_hashes" ]; then
        # Look for resources locally first
        echo "NOTICE: 1000G resources directory not found in ${HOST_DATA_DIR}/rufus_resources/kg1_hashes - downloading from S3"
        kg1_hash_arg=$($fetch_kg1_hash "$region")
    elif [ "$region" == "" ]; then
        # If we don't have a region, use entire genome wide Jhash
        kg1_hash_arg="/mnt/rufus_resources/kg1_hashes/wg.${KG1_HASH_VERSION}.Jhash"
    else
        # Convert chrN:n-m to chrN_n_m
        fmtd_reg=$(echo "$region" | tr ':-' '_')
        kg1_hash_arg="/mnt/rufus_resources/kg1_hashes/${fmtd_reg}.${KG1_HASH_VERSION}.Jhash"
    fi

    echo "$kg1_hash_arg"
}
export -f get_kg1_hash

# Looks for control hashes in ${HOST_DATA_DIR}/rufus_resources/control_hashes
# If can't find, will pull from S3
get_control_hash() {
    local region=$1
    local ctrl_hash=""

    if [ ! -d "${HOST_DATA_DIR}/rufus_resources/control_hashes" ]; then
        # Look for resources locally first
        echo "NOTICE: Internal control resources directory not found at ${HOST_DATA_DIR}/rufus_resources/control_hashes - downloading from S3"
        ctrl_hash=$($fetch_control_hash "$region")
    elif [ "$region" == "" ]; then
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

process_region() {
    local region=$1

    ctrl_arg=""
    # Check to see if controls are provided
    if [ "${#CONTROL_FILE_ARRAY[@]}" -eq 0 ]; then
        echo "No control samples provided, using internal control for single sample mode"
        internal_ctrl_hash=$($get_control_hash $region)
        ctrl_arg="-e $internal_ctrl_hash "
    else
        # Concatenate controls into a single -c delimited string
        for control in "${CONTROL_FILE_ARRAY[@]}"; do
            ctrl_arg+="-c /mnt/$control "
        done
    fi

    kg1_hash=$($get_kg1_hash $region)
    kg1_hash_arg="-e $kg1_hash"

    ref=$($get_reference)
    ref_arg="-r $ref"

    docker run -v ${HOST_DATA_DIR}:/mnt ${CONTAINER_IMAGE} bash /opt/RUFUS/runRufus.sh \
    -s /mnt/$SUBJECT_FILE \
    $ctrl_arg \
    $ref_arg \
    -m $KMER_DEPTH_CUTOFF \
    -k $KMER_LENGTH \
    -t $THREAD_LIMIT \
    $OTHER_FLAGS \
    $kg1_hash_arg \
    $region
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
echo "All RUFUS regional jobs completed. Starting merge and post-process..."

docker run -v ${HOST_DATA_DIR}:/mnt ${CONTAINER_IMAGE} bash /opt/RUFUS/post_process/post_process.sh -s "/mnt/$SUBJECT_FILE" -c "$concat_ctrl_post_arg" -r "$ref" -w "$WINDOW_SIZE" -d "/mnt"

end_time=$(date +%s)
elapsed=$((end_time - start_time))
echo "RUFUS completed. Total run time: $elapsed"