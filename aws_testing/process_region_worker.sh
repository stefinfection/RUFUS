#!/bin/bash

DEV_MOUNT="-v /home/ubuntu/RUFUS:/opt/RUFUS -v /opt/RUFUS/bin"

ENV_FILE="$1"
region="$2"

# Check RUFUS env file arg actually exists
if [ -f "$ENV_FILE" ]; then
    set -a
    source <(grep -v '^#' $ENV_FILE | grep -v '^[[:space:]]*$' | sed 's/\r$//')
    set +a
else
    echo "Error: $ENV_FILE file not found - please provide valid path to rufus.env file"
    exit 1
fi

# Check for all required variables to be filled in rufus.env
REQUIRED_VARS=(HOST_DATA_DIR SUBJECT_FILE KMER_LENGTH KMER_DEPTH_CUTOFF THREAD_LIMIT JOB_THRESHOLD REFERENCE_FASTA CONTROL_HASH_VERSION KG1_HASH_VERSION WINDOW_SIZE)
for var in "${REQUIRED_VARS[@]}"; do
    if [ -z "${!var}" ]; then
        echo "Error: Required variable $var is not set in rufus.env"
        exit 1
    fi
done

# Make sure file variables actually exist
REQUIRED_FILES=(SUBJECT_FILE REFERENCE_FASTA)
for file_var in "${REQUIRED_FILES[@]}"; do
    file_path="${HOST_DATA_DIR}/${!file_var}"
    if [ ! -f "$file_path" ]; then
        echo "Error: Required file $file_path (from variable $file_var) not found"
        exit 1
    fi
done

# Fetch functions to get resources from S3
fetch_kg1_hash() {
    local region="$1"
    local kg1_hash=""
    echo "fetching 1000G kg1 hash for region: $region" >&2

    if [ ! -d "${HOST_DATA_DIR}/rufus_resources/kg1_hashes" ]; then
        mkdir -p "${HOST_DATA_DIR}/rufus_resources/kg1_hashes"
    fi

    if [ "$region" == "" ]; then
        # If we don't have a region, use entire genome wide Jhash
        sudo docker run $DEV_MOUNT -v ${HOST_DATA_DIR}:/mnt ${RUFUS_DOCKER_IMAGE} bash /opt/RUFUS/resource_helpers/download_hash.sh "kg1" "${KG1_HASH_VERSION}" "wg"
        kg1_hash="mnt/rufus_resources/wg_kg1_${KG1_HASH_VERSION}.Jhash"
    else
        fmtd_reg=$(echo "$region" | tr ':-' '_')
        sudo docker run $DEV_MOUNT -v ${HOST_DATA_DIR}:/mnt ${RUFUS_DOCKER_IMAGE} bash /opt/RUFUS/resource_helpers/download_hash.sh "kg1" "${KG1_HASH_VERSION}" "$fmtd_reg"
        kg1_hash="/mnt/rufus_resources/${fmtd_reg}_kg1_${KG1_HASH_VERSION}.Jhash"
    fi

    echo "$kg1_hash"
}
export -f fetch_kg1_hash

fetch_control_hash() {
    local region="$1"
    local ctrl_hash=""
    echo "fetching control hash for region: $region" >&2

    if [ ! -d "${HOST_DATA_DIR}/rufus_resources/control_hashes" ]; then
        mkdir -p "${HOST_DATA_DIR}/rufus_resources/control_hashes"
    fi

    if [ "$region" == "" ]; then
        # If we don't have a region, use entire genome wide Jhash
        sudo docker run $DEV_MOUNT -v ${HOST_DATA_DIR}:/mnt ${RUFUS_DOCKER_IMAGE} bash /opt/RUFUS/resource_helpers/download_hash.sh "control" "${CONTROL_HASH_VERSION}" "wg"
        ctrl_hash="mnt/rufus_resources/wg_control_${CONTROL_HASH_VERSION}.Jhash"
    else
        # Convert chrN:n-m to chrN_n_m
        fmtd_reg=$(echo "$region" | tr ':-' '_')
        sudo docker run $DEV_MOUNT -v ${HOST_DATA_DIR}:/mnt ${RUFUS_DOCKER_IMAGE} bash /opt/RUFUS/resource_helpers/download_hash.sh "control" "${CONTROL_HASH_VERSION}" "$fmtd_reg"
        ctrl_hash="/mnt/rufus_resources/${fmtd_reg}_control_${CONTROL_HASH_VERSION}.Jhash"
    fi
    echo "$ctrl_hash"
}
export -f fetch_kg1_hash

# Checks to see if we already have BWA indexes premade for the reference argument, and points there if so
get_reference() {
    reference="${HOST_DATA_DIR}/${REFERENCE_FASTA}" 

    # If we have the exact reference BWA indexes already, point there to save some time
    if [ -d "${HOST_DATA_DIR}/rufus_resources/references" ] && [ -f "${HOST_DATA_DIR}/rufus_resources/references/${REFERENCE_FASTA}" ]; then
        reference="/mnt/rufus_resources/references/${REFERENCE_FASTA}"
    fi

    echo "$reference"
}
export -f get_reference

# Looks for 1000G hashes in ${HOST_DATA_DIR}/rufus_resources/kg1_hashes
# If can't find, will pull from S3
get_kg1_hash() {
    local region="$1"
    local kg1_hash=""

    if [ "$region" == "" ]; then
        # Use entire genome wide Jhash
        if [ -f "${HOST_DATA_DIR}/rufus_resources/kg1_hashes/wg_kg1_${KG1_HASH_VERSION}.Jhash" ]; then
            kg1_hash="/mnt/rufus_resources/kg1_hashes/wg_kg1_${KG1_HASH_VERSION}.Jhash"
        else
            echo "NOTICE: 1000G whole genome hash not found locally - downloading from S3" >&2
            kg1_hash=$(fetch_kg1_hash "$region")
        fi
    else
        # Convert chrN:n-m to chrN_n_m
        fmtd_reg=$(echo "$region" | tr ':-' '_')
        if [ -f "${HOST_DATA_DIR}/rufus_resources/kg1_hashes/${fmtd_reg}_kg1_${KG1_HASH_VERSION}.Jhash" ]; then
            kg1_hash="/mnt/rufus_resources/kg1_hashes/${fmtd_reg}_kg1_${KG1_HASH_VERSION}.Jhash"
        else
            echo "NOTICE: 1000G hash for region $region not found locally - downloading from S3" >&2
            kg1_hash=$(fetch_kg1_hash "$region")
        fi
        kg1_hash="/mnt/rufus_resources/kg1_hashes/${fmtd_reg}_kg1_${KG1_HASH_VERSION}.Jhash"
    fi

    echo "$kg1_hash"
}
export -f get_kg1_hash

# Looks for control hashes in ${HOST_DATA_DIR}/rufus_resources/control_hashes
# If can't find, will pull from S3
get_control_hash() {
    local region="$1"
    local ctrl_hash=""

    if [ "$region" == "" ]; then
        # Use entire genome wide Jhash
        if [ -f "${HOST_DATA_DIR}/rufus_resources/control_hashes/wg_control_${CONTROL_HASH_VERSION}.Jhash" ]; then
            ctrl_hash="/mnt/rufus_resources/control_hashes/wg_control_${CONTROL_HASH_VERSION}.Jhash"
        else
            echo "NOTICE: 1000G whole genome hash not found locally - downloading from S3" >&2
            ctrl_hash=$(fetch_control_hash "$region")
        fi
    else
        # Convert chrN:n-m to chrN_n_m
        fmtd_reg=$(echo "$region" | tr ':-' '_')
        if [ -f "${HOST_DATA_DIR}/rufus_resources/control_hashes/${fmtd_reg}_control_${CONTROL_HASH_VERSION}.Jhash" ]; then
            ctrl_hash="/mnt/rufus_resources/control_hashes/${fmtd_reg}_control_${CONTROL_HASH_VERSION}.Jhash"
        else
            echo "NOTICE: 1000G hash for region $region not found locally - downloading from S3" >&2
            ctrl_hash=$(fetch_control_hash "$region")
        fi
        ctrl_hash="/mnt/rufus_resources/control_hashes/${fmtd_reg}_control_${CONTROL_HASH_VERSION}.Jhash"
    fi

    echo "$ctrl_hash"
}
export -f get_control_hash


# Check to see if controls are provided
ctrl_arg=""
if [ "${#CONTROL_FILE_ARRAY[@]}" -eq 0 ]; then
    echo "No control samples provided, using internal control for single sample mode"
    internal_ctrl_hash=$(get_control_hash $region)
    ctrl_arg="-e $internal_ctrl_hash "
else
    # Concatenate controls into a single -c delimited string
    for control in "${CONTROL_FILE_ARRAY[@]}"; do
        ctrl_arg+="-c /mnt/$control "
    done
fi

kg1_hash=$(get_kg1_hash $region)
kg1_hash_arg="-e $kg1_hash"

ref=$(get_reference)

# if subject_file ends with cram, need to change region_arg to -cr
ref_arg="-r $ref"
if [[ "$SUBJECT_FILE" == *.cram ]]; then
    ref_arg="-cr $ref"
fi

region_arg=""
if [ "$region" != "" ]; then
    region_arg="-R $region"
fi


# TODO: remove mounted RUFUS code volume after testing
sudo docker run $DEV_MOUNT \
-v ${HOST_DATA_DIR}:/mnt ${RUFUS_DOCKER_IMAGE} \
bash /opt/RUFUS/runRufus.sh \
-s /mnt/$SUBJECT_FILE \
$ctrl_arg \
$ref_arg \
-m $KMER_DEPTH_CUTOFF \
-k $KMER_LENGTH \
-t $THREAD_LIMIT \
$OTHER_FLAGS \
$kg1_hash_arg \
$region_arg

# Clean up hash files
exit # TODO: remov after testing
if [ "$region" == "" ]; then
    rm ${HOST_DATA_DIR}/rufus_resources/*_hashes/wg_*.Jhash
else 
    fmtd_reg=$(echo "$region" | tr ':-' '_')
    rm ${HOST_DATA_DIR}/rufus_resources/*_hashes/$fmtd_reg*.Jhash
fi
