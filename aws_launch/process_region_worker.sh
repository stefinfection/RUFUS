#!/bin/bash

# Arguments
CONTAINER_ID="$1"
REGION="$2"

# Constants
ENV_FILE="/mnt/rufus_resources/rufus.env"
DEFAULT_KG1_HASH_VERSION="v3.0"
DEFAULT_CONTROL_HASH_VERSION="v1.0"

# Import env file now that we're inside of container
set -a
source <(grep -v '^#' $ENV_FILE | grep -v '^[[:space:]]*$' | sed 's/\r$//')
set +a

# Fetches control or kg1 hash from S3 for region if region arg provided, or whole genome hash otherwise
# Returns path to downloaded hash
# Fetches control or kg1 hash from S3 for region if region arg provided, or whole genome hash otherwise
# Returns path to downloaded hash
fetch_hash() {
    local region="$1"
    local hash_type="$2"
    # Capitalize all letters in hash_type
    local hash_type_upper=$(echo "$hash_type" | tr '[:lower:]' '[:upper:]')
    # Output var
    local hash=""
    
    # Get the hash version - use specific version if set, otherwise use default
    local version_var="${hash_type_upper}_HASH_VERSION"
    local hash_version="${!version_var}"
    
    if [ -z "$hash_version" ]; then
        local default_var="DEFAULT_${hash_type_upper}_HASH_VERSION"
        hash_version="${!default_var}"
    fi
    
    if [ -z "$region" ]; then
        # If we don't have a region, use entire genome wide Jhash
        echo "Fetching version ${hash_version} whole genome ${hash_type} hash" >&2
        docker exec ${CONTAINER_ID} bash /opt/RUFUS/resource_helpers/download_hash.sh "${hash_type}" "${hash_version}" "wg"
        hash="/mnt/rufus_resources/${hash_type}_hashes/wg_${hash_type}_${hash_version}.Jhash"
    else
        # Convert chrN:n-m to chrN_n_m
        echo "Fetching version ${hash_version} ${hash_type} hash for region: $region" >&2
        local fmtd_reg=$(echo "$region" | tr ':-' '_')
        docker exec ${CONTAINER_ID} bash /opt/RUFUS/resource_helpers/download_hash.sh "${hash_type}" "${hash_version}" "$fmtd_reg"
        hash="/mnt/rufus_resources/${hash_type}_hashes/${fmtd_reg}_${hash_type}_${hash_version}.Jhash"
    fi
    
    echo "$hash"
}
export -f fetch_hash

# Looks for control or kg1 hashes in /mnt/rufus_resources/${hash_type}_hashes first
# At this point, we know that if a local directory has been provided and mounted, it contains at least one *.Jhash file
# If local directory contains multiple *.Jhash files matching region or WG, will return error code
# If can't find, will pull from S3
get_hash() {
    local region="$1"
    local geo_type="$2"
    local hash_type="$3"

    local hash_type_upper=$(echo "$hash_type" | tr '[:lower:]' '[:upper:]')

    # Output vars
    local hash=""
    local output_echo=""
    
    # Get the user's original directory path for error messages
    local env_var="${hash_type_upper}_HASH_LOCAL_DIR"
    local user_dir="${!env_var}"
    
    # The actual mounted directory in the container
    local local_dir="/mnt/rufus_resources/${hash_type}_hashes"
    
    # Local hashes
    if [ "$geo_type" == "local" ]; then
        # Whole genome mode and we're looking locally
        if [ "$region" == "" ]; then
            file_count=$(ls "${local_dir}"/*wg*${hash_type}*.Jhash 2>/dev/null | wc -l)
            if [ "$file_count" -gt 1 ]; then
                echo "ERROR: Multiple whole genome ${hash_type} hash files found in ${user_dir}:" >&2
                ls "${local_dir}"/*wg*${hash_type}*.Jhash >&2
                echo "Please ensure only one *wg*${hash_type}*.Jhash file exists in the directory." >&2
                return 1
            elif [ "$file_count" -eq 1 ]; then
                hash=$(ls "${local_dir}"/*wg*${hash_type}*.Jhash)
                output_echo+="Using local whole genome ${hash_type} hash at $hash"
            else
                # File does not exist locally, fetch from S3
                output_echo+="Could not find local whole genome ${hash_type} hash in ${user_dir}. The file in this directory must be named like \"*wg*${hash_type}*.Jhash\" for RUFUS to recognize it. Will attempt to fetch from S3."
                hash=$(fetch_hash "$region" "$hash_type")
            fi
        # Region mode and we're looking locally
        else
            fmtd_reg=$(echo "$region" | tr ':-' '_')
            file_count=$(ls "${local_dir}"/*${fmtd_reg}*${hash_type}*.Jhash 2>/dev/null | wc -l)
            if [ "$file_count" -gt 1 ]; then
                echo "ERROR: Multiple ${hash_type} hash files for region $region found in ${user_dir}:" >&2
                ls "${local_dir}"/*${fmtd_reg}*${hash_type}*.Jhash >&2
                echo "Please ensure only one *${fmtd_reg}*${hash_type}*.Jhash file exists in the directory." >&2
                return 1
            elif [ "$file_count" -eq 1 ]; then
                hash=$(ls "${local_dir}"/*${fmtd_reg}*${hash_type}*.Jhash)
                output_echo+="Using local ${hash_type} hash for region $region at: $hash"
            else
                # File does not exist locally, fetch from S3
                output_echo+="Could not find local ${hash_type} hash for region $region in ${user_dir}. The file in this directory must be named like \"*${fmtd_reg}*${hash_type}*.Jhash\" for RUFUS to recognize it. Will attempt to fetch from S3."
                hash=$(fetch_hash "$region" "$hash_type")
            fi
        fi
    # Remote fetching of hashes
    else
        output_echo+="Fetching ${hash_type} hash from S3."
        hash=$(fetch_hash "$region" "$hash_type")
    fi
    
    echo -e "$output_echo" >&2
    echo "$hash"
}
export -f get_hash


# Compose control argument of both or one of control hashes and paired control files
ctrl_arg=""
if [ "$CONTROL_HASH_LOCAL_DIR" != "" ]; then
    control_hash=$(get_hash $REGION "local" "control") || exit 1
    ctrl_arg+="-e $control_hash "
elif [ "$CONTROL_HASH_VERSION" != "" ]; then
    control_hash=$(get_hash $REGION "remote" "control") || exit 1
    ctrl_arg+="-e $control_hash "
elif [ "${#CONTROL_FILE_ARRAY[@]}" -eq 0 ]; then
    echo "No local control hashes, paired controls, or control hash version provided, fetching default $DEFAULT_CONTROL_HASH_VERSION hashes piecemeal" >&2
    CONTROL_HASH_VERSION="$DEFAULT_CONTROL_HASH_VERSION"
    control_hash=$(get_hash $REGION "remote" "control") || exit 1
    ctrl_arg+="-e $control_hash "
fi

# If we have paired controls provided, also use those
if [ "${#CONTROL_FILE_ARRAY[@]}" -ne 0 ]; then
    # Concatenate controls into -c delimited string
    for control in "${CONTROL_FILE_ARRAY[@]}"; do
        control_base=$(basename "$control")
        ctrl_arg+="-c /mnt/$control_base "
    done
fi

# Compose kg1 hash argument
if [ "$KG1_HASH_LOCAL_DIR" != "" ]; then
    kg1_hash=$(get_hash $REGION "local" "kg1") || exit 1
    kg1_hash_arg="-e $kg1_hash"
elif [ "$KG1_HASH_VERSION" != "" ]; then
    kg1_hash=$(get_hash $REGION "remote" "kg1") || exit 1
    kg1_hash_arg="-e $kg1_hash"
fi

ref_base=$(basename "$REFERENCE_FASTA")
ref_arg="-r /mnt/$ref_base"
# If subject_file ends with cram, need to change region_arg to -cr
subject_base=$(basename "$SUBJECT_FILE")
if [[ "$subject_base" == *.cram ]]; then
    ref_arg="-cr /mnt/$ref_base"
fi

# Region arg
if [ "$REGION" != "" ]; then
    region_arg="-R $REGION"
    fmtd_reg=$(echo "$REGION" | tr ':-' '_')
fi

RUFUS_CMD="/opt/RUFUS/runRufus.sh \
  -s /mnt/$subject_base \
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
if [ "$REGION" == "" ]; then
    if [ "$KEEP_DOWNLOADED_HASHES" != "TRUE" ] && [ "$KEEP_DOWNLOADED_HASHES" != "true" ]; then
        docker exec ${CONTAINER_ID} rm /mnt/rufus_resources/control_hashes/*wg*control*.Jhash
        docker exec ${CONTAINER_ID} rm /mnt/rufus_resources/kg1_hashes/*wg*kg1*.Jhash
    fi
else 
    fmtd_reg=$(echo "$REGION" | tr ':-' '_')
    if [ "$KEEP_DOWNLOADED_HASHES" != "TRUE" ] && [ "$KEEP_DOWNLOADED_HASHES" != "true" ]; then
        docker exec ${CONTAINER_ID} rm /mnt/rufus_resources/control_hashes/*$fmtd_reg*.Jhash
        docker exec ${CONTAINER_ID} rm /mnt/rufus_resources/kg1_hashes/*$fmtd_reg*.Jhash
    fi
fi
