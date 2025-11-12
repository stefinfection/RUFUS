#!/bin/bash

ENV_FILE="$1"
CONTAINER_ID="$2"
region="$3"

# Constants
DEFAULT_KG1_HASH_VERSION="v3.0"
DEFAULT_CONTROL_HASH_VERSION="v1.0"

# Globals
delete_kg1_hashes=false
delete_control_hashes=false


# Check RUFUS env file arg actually exists
if [ -f "$ENV_FILE" ]; then
    set -a
    source <(grep -v '^#' $ENV_FILE | grep -v '^[[:space:]]*$' | sed 's/\r$//')
    set +a
else
    echo "Error: $ENV_FILE file not found - please provide valid path to rufus.env file"
    exit 1
fi



# Check that we have at least one of the following filled out: KG1_HASH_VERSION KG1_HASH_LOCAL_DIR
if [ -z "$KG1_HASH_VERSION" ] && [ -z "$KG1_HASH_LOCAL_DIR" ]; then
    echo "Error: Must provide at least one of: KG1_HASH_VERSION or KG1_HASH_LOCAL_DIR in rufus.env"
    exit 1
else
    if [ "$KG1_HASH_LOCAL_DIR" != "" ]; then
        echo "Using local 1000G hashes at: $KG1_HASH_LOCAL_DIR"
    elif [ "$KG1_HASH_VERSION" != "" ]; then
        echo "Fetching S3 1000G hashes with version: $KG1_HASH_VERSION"
    fi
fi

# Fetch functions to get resources from S3
# TODO: add in option to pass flag here to keep hashes?
fetch_kg1_hash() {
    local region="$1"
    local kg1_hash=""
    echo "fetching 1000G kg1 hash for region: $region" >&2
    delete_kg1_hashes=true

    # We may not have a hash version if we're missing one of the hashes in our local dir and don't have KG1_HASH_VERSION set
    local hash_version="$KG1_HASH_VERSION"
    if [ "$hash_version" == "" ]; then
        # Set to the latest default for now
        hash_version="v3.0"
    fi

    hash_dir=""
    if [ "${KG1_HASH_LOCAL_DIR}" == "" ]; then
        hash_dir="${WORKING_DIR}/rufus_resources/kg1_hashes"
        if [ ! -d "$hash_dir" ]; then
            mkdir -p "$hash_dir"
        fi
    else
        hash_dir="${KG1_HASH_LOCAL_DIR}"
    fi

    if [ "$region" == "" ]; then
        # If we don't have a region, use entire genome wide Jhash
        docker exec ${CONTAINER_ID} bash /opt/RUFUS/resource_helpers/download_hash.sh "kg1" "${hash_version}" "wg"
        kg1_hash="mnt/rufus_resources/wg_kg1_${hash_version}.Jhash"
    else
        fmtd_reg=$(echo "$region" | tr ':-' '_')
        docker exec ${CONTAINER_ID} bash /opt/RUFUS/resource_helpers/download_hash.sh "kg1" "${KG1_HASHhash_version_VERSION}" "$fmtd_reg"
        kg1_hash="/mnt/rufus_resources/${fmtd_reg}_kg1_${hash_version}.Jhash"
    fi

    echo "$kg1_hash"
}
export -f fetch_kg1_hash

# TODO: add in option to pass flag here to keep hashes?
fetch_control_hash() {
    local region="$1"
    local ctrl_hash=""
    echo "fetching control hash for region: $region" >&2
    delete_control_hashes=true

    # We may not have a hash version if we're missing one of the hashes in our local dir and don't have CONTROL_HASH_VERSION set
    local hash_version="$CONTROL_HASH_VERSION"
    if [ "$hash_version" == "" ]; then
        # Set to the latest default for now
        hash_version="v1.0"
    fi

    if [ ! -d "${WORKING_DIR}/rufus_resources/control_hashes" ]; then
        mkdir -p "${WORKING_DIR}/rufus_resources/control_hashes"
    fi

    if [ "$region" == "" ]; then
        # If we don't have a region, use entire genome wide Jhash
        docker exec ${CONTAINER_ID} bash /opt/RUFUS/resource_helpers/download_hash.sh "control" "${hash_version}" "wg"
        ctrl_hash="mnt/rufus_resources/wg_control_${hash_version}.Jhash"
    else
        # Convert chrN:n-m to chrN_n_m
        fmtd_reg=$(echo "$region" | tr ':-' '_')
        docker exec ${CONTAINER_ID} bash /opt/RUFUS/resource_helpers/download_hash.sh "control" "${hash_version}" "$fmtd_reg"
        ctrl_hash="/mnt/rufus_resources/${fmtd_reg}_control_${hash_version}.Jhash"
    fi
    echo "$ctrl_hash"
}
export -f fetch_control_hash

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

# Looks for 1000G hashes in ${KG1_HASH_LOCAL_DIR}
# If can't find, will pull from S3
get_kg1_hash() {
    local region="$1"
    local kg1_hash=""
    local kg1_output_clause=""

    if [ "$KG1_HASH_LOCAL_DIR" != "" ]; then

        if [ "$region" == "" ]; then
            # We don't have a region and we're looking locally
            if [ -d "${KG1_HASH_LOCAL_DIR}" ] && ls "$KG1_HASH_LOCAL_DIR"/*wg*kg1*.Jhash 1> /dev/null 2>&1; then
                # File exists
                kg1_output_clause+="Using local whole genome 1000G hash at: $KG1_HASH_LOCAL_DIR"
                kg1_base=$(ls "$KG1_HASH_LOCAL_DIR"/*wg*kg1*.Jhash)
                kg1_hash="/mnt/rufus_resources/$kg1_base"
            else
                # File does not exist locally, fetch from S3
                kg1_output_clause+=" Could not find local whole genome 1000G hash at: $KG1_HASH_LOCAL_DIR, will attempt to fetch from S3"
                kg1_hash=$(fetch_kg1_hash "$region") 
            fi
        else
            # We have a region and we're looking locally
            fmtd_reg=$(echo "$region" | tr ':-' '_')
            # Look for a file with our formatted region in the name (because we don't enforce version if looking locally)
            if ls "$KG1_HASH_LOCAL_DIR"/*${fmtd_reg}*kg1*.Jhash 1> /dev/null 2>&1; then
                kg1_output_clause+="Using local 1000G hash for region $region at: $KG1_HASH_LOCAL_DIR"
                kg1_file=$(ls "$KG1_HASH_LOCAL_DIR"/*${fmtd_reg}*kg1*.Jhash)
                kg1_base=$(basename "$kg1_file")
                kg1_hash="/mnt/rufus_resources/$kg1_base"
            else
                # File does not exist locally, fetch from S3
                kg1_output_clause+=" Could not find local 1000G hash for region $region at: $KG1_HASH_LOCAL_DIR, will attempt to fetch from S3"
                kg1_hash=$(fetch_kg1_hash "$region") 
            fi
        fi
    elif [ "$KG1_HASH_VERSION" != "" ]; then
    # Then use S3 if we don't have local hashes
        if [ "$region" == "" ]; then
            kg1_output_clause+="Fetching S3 whole genome 1000G hash with version: $KG1_HASH_VERSION"
            kg1_hash=$(fetch_kg1_hash "$region")
        else
            kg1_output_clause+="Fetching S3 1000G hash for region $region with version: $KG1_HASH_VERSION"
            kg1_hash=$(fetch_kg1_hash "$region")
        fi
    fi
    echo $kg1_output_clause >&2
    echo "$kg1_hash"
}
export -f get_kg1_hash

# Looks for control hashes in user provided ${CONTROL_HASH_LOCAL_DIR}
# If can't find, will pull from S3
get_control_hash() {
    local region="$1"
    local type="$2"

    # Output vars
    local ctrl_hash=""
    local dir_mount=""
    local ctrl_output_echo=""

    # Local hashes
    if [ type == "local" ]; then

        # Whole genome mode and we're looking locally
        if [ "$region" == "" ]; then
            curr_hash=$(ls "$CONTROL_HASH_LOCAL_DIR"/*wg*control*.Jhash)
            if [ "$curr_hash" != "" ]; then
                control_output_echo+="Using local whole genome control hash at: $curr_hash"
                curr_base=$(basename "$curr_hash")
                ctrl_hash="/mnt/rufus_resources/$curr_base"
                dir_mount="${CONTROL_HASH_LOCAL_DIR}:/mnt/rufus_resources/control_hashes"
            else
                # File does not exist locally, fetch from S3
                control_output_echo+=" Could not find local whole genome control hash in: $CONTROL_HASH_LOCAL_DIR. This file must be formatted like "*wg*control*.Jhash" for RUFUS to recognize it. Will attempt to fetch from S3."
                ctrl_hash=$(fetch_control_hash "$region")
            fi
        # Region mode and we're looking locally
        else
            fmtd_reg=$(echo "$region" | tr ':-' '_')
            curr_hash=$(ls "$CONTROL_HASH_LOCAL_DIR"/*${fmtd_reg}*control*.Jhash)
            if [ "$curr_hash" != "" ]; then
                control_output_echo+="Using local control hash for region $region at: $curr_hash"
                curr_base=$(basename "$curr_hash")
                ctrl_hash="/mnt/rufus_resources/$curr_base"
            else
                # File does not exist locally, fetch from S3
                control_output_echo+=" Could not find local control hash for region $region in: $CONTROL_HASH_LOCAL_DIR. This file must be formatted like "chrN_start_end*control*.Jhash" for RUFUS to recognize it. Will attempt to fetch from S3."
                ctrl_hash=$(fetch_control_hash "$region") 
            fi
        fi
    # Remote fetching of hashes
    else
        if [ "$region" == "" ]; then
            control_output_echo+="Fetching S3 whole genome control hash with version: $CONTROL_HASH_VERSION."
            ctrl_hash=$(fetch_control_hash "$region")
        else
            control_output_echo+="Fetching S3 control hash for region $region with version: $CONTROL_HASH_VERSION."
            ctrl_hash=$(fetch_control_hash "$region")
        fi
    fi
    echo $control_output_echo >&2
    echo "$ctrl_hash"
}
export -f get_control_hash


# Compose control argument of both or one of control hashes and paired control files
ctrl_arg=""

if [ "$CONTROL_HASH_LOCAL_DIR" != "" ]; then
    control_hash=$(get_control_hash $region "local")
    ctrl_arg+="-e $control_hash "
elif [ "$CONTROL_HASH_VERSION" != "" ]; then
    control_hash=$(get_control_hash $region "remote")
    ctrl_arg+="-e $control_hash "
elif [ "${#CONTROL_FILE_ARRAY[@]}" -eq 0 ]; then
    echo "No control hash version provided, using default: $DEFAULT_CONTROL_HASH_VERSION" >&2
    CONTROL_HASH_VERSION="$DEFAULT_CONTROL_HASH_VERSION"
    control_hash=$(get_control_hash $region "remote")
    ctrl_arg+="-e $control_hash "
fi

# If we have paired controls provided, also use those
if [ "${#CONTROL_FILE_ARRAY[@]}" -ne 0 ]; then
    # Concatenate controls into -c delimited string
    for control in "${CONTROL_FILE_ARRAY[@]}"; do
        ctrl_arg+="-c /mnt/$control "
    done
fi

if [ "$KG1_HASH_VERSION" == "" ] && [ "$KG1_HASH_LOCAL_DIR" == "" ]; then
    echo "WARNING: not removing population variants found in the 1000G cohort as both KG1_HASH_VERSION and KG1_HASH_LOCAL_DIR are empty" >&2
else
    kg1_hash=$(get_kg1_hash $region)
    kg1_hash_arg="-e $kg1_hash"
fi

# If subject_file ends with cram, need to change region_arg to -cr
ref=$(get_reference)
ref_arg="-r $ref"
subject_base=$(basename "$SUBJECT_FILE")
if [[ "$subject_base" == *.cram ]]; then
    ref_arg="-cr $ref"
fi

# Region args
region_arg=""
fmtd_reg=""
if [ "$region" != "" ]; then
    region_arg="-R $region"
    fmtd_reg=$(echo "$region" | tr ':-' '_')
fi

# Make log directory
mkdir -p ${HOST_DATA_DIR}/rufus_resources/logs

# TODO: will need to mount any input file directories
RUFUS_CMD="/opt/RUFUS/runRufus.sh \
  -s /mnt/$SUBJECT_FILE \
  $ctrl_arg \
  $ref_arg \
  -m $KMER_DEPTH_CUTOFF \
  -k $KMER_LENGTH \
  -t $THREAD_LIMIT \
  $OTHER_FLAGS \
  $kg1_hash_arg \
  $region_arg"

docker exec "$CONTAINER_ID" bash -c \
    "$RUFUS_CMD \
    > /mnt/rufus_resources/logs/${fmtd_reg}.out \
    2> /mnt/rufus_resources/logs/${fmtd_reg}.err"

# Clean up hash files
if [ "$region" == "" ]; then
    if $delete_control_hashes; then
        rm ${HOST_DATA_DIR}/rufus_resources/control_hashes/wg_control*.Jhash
    fi

    if $delete_kg1_hashes; then
        rm ${HOST_DATA_DIR}/rufus_resources/kg1_hashes/wg_kg1*.Jhash
    fi
else 
    fmtd_reg=$(echo "$region" | tr ':-' '_')

    if $delete_control_hashes; then
        rm ${HOST_DATA_DIR}/rufus_resources/control_hashes/$fmtd_reg*.Jhash
    fi

    if $delete_kg1_hashes; then
        rm ${HOST_DATA_DIR}/rufus_resources/kg1_hashes/$fmtd_reg*.Jhash
    fi
fi
