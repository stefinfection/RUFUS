#!/bin/bash

# Fill in
HOST_DATA_DIR=


start_time=$(date +%s)

# Check for RUFUS env file
if [ -f "$HOST_DATA_DIR/rufus_resources/rufus.env" ]; then
    set -a
    source <(grep -v '^#' $PATH_TO_ENV/rufus.env | grep -v '^[[:space:]]*$' | sed 's/\r$//')
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

get_reference() {
    if [ ! -d "${HOST_DATA_DIR}/rufus_resources/references" ] || [ ! -f "${HOST_DATA_DIR}/rufus_resources/references/${REFERENCE_FASTA}" ]; then
        # Check for resources directory and that the reference is in there - warn if not
        echo -n "Warning - ${REFERENCE_FASTA} not found in ${HOST_DATA_DIR}/rufus_resources/references directory - "
        echo "Utilizing pre-built BWA references provided in the RUFUS resources directory can speed up run time next time"
        reference="mnt/${REFERENCE_FASTA}"
    else
        # If it is in references, use it 
        reference="/mnt/rufus_resources/references/${REFERENCE_FASTA}"
    fi

    echo "$reference"
}

# Fetch resources from S3 if not present locally
fetch_kg1_hash() {
    local region=$1

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
        fmtd_reg=$(echo "$region" | tr ':-' '_')
        chrom=$(echo "$fmtd_reg" | cut -d'_' -f1)
        aws s3 sync "s3://rufus.marth.lab/public_access_data/kg1_hashes/${KG1_HASH_VERSION}/${chrom}/${fmtd_reg}_kg1_${KG1_HASH_VERSION}.Jhash" "${HOST_DATA_DIR}/rufus_resources/${fmtd_reg}_kg1_${KG1_HASH_VERSION}.Jhash"
        kg1_hash="/mnt/rufus_resources/${fmtd_reg}_kg1_${KG1_HASH_VERSION}.Jhash"
    fi

    echo "$kg1_hash"
}

fetch_control_hash() {
    local region=$1

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
        fmtd_reg=$(echo "$region" | tr ':-' '_')
        chrom=$(echo "$fmtd_reg" | cut -d'_' -f1)
        aws s3 sync "s3://rufus.marth.lab/public_access_data/control_hashes/${CONTROL_HASH_VERSION}/${chrom}/${fmtd_reg}_control_${CONTROL_HASH_VERSION}.Jhash" "${HOST_DATA_DIR}/rufus_resources/${fmtd_reg}_control_${CONTROL_HASH_VERSION}.Jhash"
        ctrl_hash="/mnt/rufus_resources/${fmtd_reg}_control_${CONTROL_HASH_VERSION}.Jhash"
    fi
    echo "$ctrl_hash"
}

process_region() {
    local region=$1

    ctrl_arg=""
    # Check to see if controls are provided
    if [ "${#CONTROL_FILE_ARRAY[@]}" -eq 0 ]; then
        echo "No control samples provided, using internal control for single sample mode"
        internal_ctrl_hash=$($fetch_control_hash $region)
        ctrl_arg="-e $internal_ctrl_hash "
    else
        # Concatenate controls into a single -c delimited string
        for control in "${CONTROL_FILE_ARRAY[@]}"; do
            ctrl_arg+="-c /mnt/$control "
        done
    fi

    kg1_hash=$(fetch_kg1_hash $region)
    kg1_hash_arg="-e $kg1_hash"

    ref=$(get_reference)
    ref_arg="-r $ref"

    docker run --rm -v ${HOST_DATA_DIR}:/mnt ${CONTAINER_IMAGE} bash /opt/RUFUS/runRufus.sh \
    -s /mnt/$SUBJECT_FILE \
    $ctrl_arg \
    $ref_arg \
    -m $KMER_DEPTH_CUTOFF \
    -k $KMER_LENGTH \
    -t $THREAD_LIMIT \
    $OTHER_FLAGS \
    $kg1_hash_arg \
    $region

    # Uncomment below to use singularity
    # singularity exec --bind ${HOST_DATA_DIR}:/mnt ${CONTAINER_IMAGE} bash /opt/RUFUS/runRufus.sh \
    #     -s /mnt/$SUBJECT_FILE \
    #     $ctrl_arg \
    #     $ref_arg \
    #     -m $KMER_DEPTH_CUTOFF \
    #     -k $KMER_LENGTH \
    #     -t $THREAD_LIMIT \
    #     $OTHER_FLAGS \
    #     $kg1_hash_arg \
    #     $region
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

docker run --rm -v ${HOST_DATA_DIR}:/mnt ${CONTAINER_IMAGE} bash /opt/RUFUS/post_process/post_process.sh -s "/mnt/$SUBJECT_FILE" -c "$concat_ctrl_post_arg" -r "$ref" -w "$WINDOW_SIZE" -d "/mnt"

# Uncomment below to use singularity
#singularity exec --bind "${HOST_DATA_DIR}:/mnt" "${CONTAINER_IMAGE}$" bash /opt/RUFUS/post_process/post_process.sh -s "/mnt/$SUBJECT_FILE" -c "$concat_ctrl_post_arg" -r "$ref" -w "$WINDOW_SIZE" -d "/mnt"

end_time=$(date +%s)
elapsed=$((end_time - start_time))
echo "RUFUS completed. Total run time: $elapsed"
