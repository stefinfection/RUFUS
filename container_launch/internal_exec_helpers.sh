#!/bin/bash
# A collection of functions run inside the container to get command line arguments 
# for RUFUS and post-process execution.
# Requires the RUFUS config yaml to be parsed and translated into a cleaned env file.
# SJ Georges Nov2025

# ----------------- EXEC FUNCTIONS ----------------- #

# TODO: do I need to resource env file here?

# TODO: left off here - need to adapt all of these functions like mount section
# should I return string array for these and mount_clause below and then printf in write functions?

# Fetches control or kg1 hash from S3 for region if region arg provided, or whole genome hash otherwise
# Returns path inside container to downloaded hash directory (named REMOTE_CONTROL/KG1_HASH_DIR in globals.env)
# Hash type
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

    local remote_dir="INTERNAL_REMOTE_${hash_type_upper}_HASH_DIR"
    local remote_hash_dir="${!remote_dir}"
    
    if [ -z "$hash_version" ]; then
        local default_var="DEFAULT_${hash_type_upper}_HASH_VERSION"
        hash_version="${!default_var}"
    fi
    
    if [ -z "$region" ]; then
        # If we don't have a region, use entire genome wide Jhash
        echo "Fetching version ${hash_version} whole genome ${hash_type} hash" >&2
        docker exec ${CONTAINER_ID} bash "$DOWNLOAD_HASH_SCRIPT" "${hash_type}" "${hash_version}" "wg" >&2
        hash="${remote_hash_dir}wg_${hash_type}_${hash_version}.Jhash"
    else
        # Convert chrN:n-m to chrN_n_m
        echo "Fetching version ${hash_version} ${hash_type} hash for region: $region" >&2
        local fmtd_reg=$(echo "$region" | tr ':-' '_')
        docker exec ${CONTAINER_ID} bash "${DOWNLOAD_HASH_SCRIPT}" "${hash_type}" "${hash_version}" "$fmtd_reg" >&2
        hash="${remote_hash_dir}${fmtd_reg}_${hash_type}_${hash_version}.Jhash"
    fi

    echo "$hash"
}
export -f fetch_hash

# Looks for control or kg1 hashes (in $LOCAL_{CONTROL/KG1}_HASH_DIR in globals.env) first
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
    local host_dir="${!env_var}"
    local cont_dir="${RUNTIME_TEMP_DIR}${hash_type}_hashes"

    # Local hashes
    if [ "$geo_type" == "local" ]; then
        # Whole genome mode and we're looking locally
        if [ "$region" == "" ]; then
            file_count=$(docker exec "$CONTAINER_ID" bash -c "ls ${cont_dir}/*wg*.Jhash 2>/dev/null | wc -l")
            if [ "$file_count" -gt 1 ]; then
                echo "ERROR: Multiple whole genome ${hash_type} hash files found in ${host_dir}:" >&2
                ls "${host_dir}"/*wg*.Jhash >&2
                echo "Please ensure only one *wg*.Jhash file exists in the directory." >&2
                return 1
            elif [ "$file_count" -eq 1 ]; then
                hash=$(docker exec "$CONTAINER_ID" bash -c "ls ${cont_dir}/*wg*.Jhash")
                output_echo+="Using local whole genome ${hash_type} hash at $hash"
            else
                # File does not exist locally, fetch from S3
                output_echo+="Could not find local whole genome ${hash_type} hash in ${host_dir}. The file in this directory must be named like \"*wg*.Jhash\" for RUFUS to recognize it. Will attempt to fetch from S3."
                hash=$(fetch_hash "$region" "$hash_type")
            fi
        # Region mode and we're looking locally
        else
            fmtd_reg=$(echo "$region" | tr ':-' '_')
            file_count=$(docker exec "$CONTAINER_ID" bash -c "ls ${cont_dir}/*${fmtd_reg}*.Jhash 2>/dev/null | wc -l")
            if [ "$file_count" -gt 1 ]; then
                echo "ERROR: Multiple ${hash_type} hash files for region $region found in ${host_dir}:" >&2
                ls "${host_dir}"/*${fmtd_reg}*.Jhash >&2
                echo "Please ensure only one *${fmtd_reg}*.Jhash file exists in the directory." >&2
                return 1
            elif [ "$file_count" -eq 1 ]; then
                hash=$(docker exec "$CONTAINER_ID" bash -c "ls ${cont_dir}/*${fmtd_reg}*.Jhash")
                output_echo+="Using local ${hash_type} hash for region $region at: $hash"
            else
                # File does not exist locally, fetch from S3
                output_echo+="Could not find local ${hash_type} hash for region $region in ${host_dir}. The file in this directory must be named like \"*${fmtd_reg}*.Jhash\" for RUFUS to recognize it. Will attempt to fetch from S3."
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

# Returns RUFUS argument string for paired controls and prebuilt hashes as appropriate
# Assumes required args met and parsing into temp env file completed
get_rufus_control_arg() {
    ctrl_arg=""
    if [ "$CONTROL_HASH_LOCAL_DIR" != "" ]; then
        control_hash=$(get_hash $REGION "local" "control") || exit 1
        ctrl_arg+="-e $control_hash "
    elif [ "$CONTROL_HASH_VERSION" != "" ]; then
        echo "Using control hash version: $CONTROL_HASH_VERSION" >&2
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
            ctrl_arg+="-c ${RUNTIME_TEMP_DIR}${control_base} "
        done
    fi

}
export -f get_rufus_control_arg

# TODO: only call this if override skip arg not set
# Returns RUFUS argument string for 1000G prebuilt hashes as appropriate
# Assumes required args met and parsing into temp env file completed
get_rufus_kg1_arg() {
    if [ "$KG1_HASH_LOCAL_DIR" != "" ]; then
        kg1_hash=$(get_hash $REGION "local" "kg1") || exit 1
        kg1_hash_arg="-e $kg1_hash"
    elif [ "$KG1_HASH_VERSION" != "" ]; then
        kg1_hash=$(get_hash $REGION "remote" "kg1") || exit 1
        kg1_hash_arg="-e $kg1_hash"
    fi

    echo "$kg1_hash"
}
export -f get_rufus_kg1_arg

# Returns RUFUS argument string for reference
# Agnostic to BWA index status
get_rufus_ref_arg() {
    ref_base=$(basename "$REFERENCE_FASTA")
    ref_arg="-r ${REF_INDEX_DIR}${ref_base}"

    # If subject_file ends with cram, need to change region_arg to -cr
    subject_base=$(basename "$SUBJECT_FILE")
    if [[ "$subject_base" == *.cram ]]; then
        ref_arg="-cr ${REF_INDEX_DIR}${ref_base}"
    fi
}
export -f get_rufus_ref_arg

# Returns single string of arguments provided directly to RUFUS run script
# All args here are not relative to a region
# TODO: do I want to return array of strings here and do printf at PoC to ensure proper formatting?
get_rufus_args() {
    local rufus_args=""

    subject_base=$(basename "${SUBJECT_FILE}")
    ref_base=$(basename "${REFERENCE_FASTA}")
    ref_arg="-r ${REF_INDEX_DIR}$ref_base"

    if [ "${#CONTROL_FILE_ARRAY[@]}" -ne 0 ]; then
        # Concatenate controls into -c delimited string
        for control in "${CONTROL_FILE_ARRAY[@]}"; do
            control_base=$(basename "$control")
            ctrl_arg+="-c ${RUNTIME_TEMP_DIR}${control_base} "
        done
    fi

    rufus_args+="\
    -s ${RUNTIME_TEMP_DIR}${subject_base} \
    $ctrl_arg \
    $ref_arg \
    -m $KMER_DEPTH_CUTOFF \
    -k $KMER_LENGTH \
    -t $THREAD_LIMIT \
    $OTHER_FLAGS"

    echo "$rufus_args"
}
export -f get_rufus_args

get_post_process_args() {
    subject_base=$(basename "${SUBJECT_FILE}")
    ref_base=$(basename "${REFERENCE_FASTA}")

    concat_ctrl_post_arg=""
    if [ ${#CONTROL_FILE_ARRAY[@]} -gt 0 ]; then
        # Concatenate controls without -c delimiters
        concat_ctrls=""
        for control in "${CONTROL_FILE_ARRAY[@]}"; do
            ctrl_base=$(basename "${control}")
            concat_ctrls+="${RUNTIME_TEMP_DIR}${ctrl_base} "
        done
        concat_ctrl_post_arg="-c $concat_ctrls"
    fi

    echo "-s ${RUNTIME_TEMP_DIR}${subject_base} \
        -r ${REF_INDEX_DIR}${ref_base} \
        -w $WINDOW_SIZE \
        -d $RUNTIME_TEMP_DIR \
        $concat_ctrl_post_arg"
}
export -f get_post_process_args

# ----------------- OTHER HELPERS ----------------- #

# Takes in /temp/rufus_config.yaml
# Parses and error checks arguments
# Writes arguments to temp.env
parse_config_file() {
    # TODO: when create cleaned env file, also source globals so we have access
    # Globals will have $CONFIG_PATH that parser should have access to
    # call config parser.py - do all parsing and input checking here
    # have it write to temp_env file in safe way
}

# Returns TRUE if we don't have all of the BWA indexes of our reference file in the same directory as the original
check_for_bwa_indexes() {
    local build_refs="FALSE"
    local ref_base=""
    ref_base=$(basename "$REFERENCE_FASTA")
    local internal_ref_path="/mnt/rufus_temp/bwa_indexes/$ref_base"

    # Determine the base filename (without .gz if present)
    if [[ "$internal_ref_path" == *.gz ]]; then
        ref_file="${internal_ref_path%.gz}"
    else
        ref_file="$internal_ref_path"
    fi

    for ext in sa bwt pac amb ann fai; do
        if [[ ! -e "${ref_file}.${ext}" ]]; then
            build_refs="TRUE"
            break
        fi
    done

    echo "$build_refs"
}
export -f check_for_bwa_indexes

# Returns slurm specs for regional job
get_slurm_specs() {

    # TODO: how to consolidate if user provides thread limit that doesn't make sense with cpus-per-task

    # Default is split jobs evenly on num nodes provided
    local num_regions=$(wc -l "$REGION_FILE")
    local ntasks=$($num_regions/$SLURM_NODES) # TODO: need to account for uneven job #

    local cpus_per_task="10"
    if [ -n "$SLURM_CPUS_PER_TASK" ]; then
        cpus_per_task=$SLURM_CPUS_PER_TASK
    fi

    local mem_per_task="8G"
    if [ -n "$SLURM_MEM_PER_TASK" ]; then
        mem_per_task=$SLURM_MEM_PER_TASK
    fi

    echo "$ntasks $cpus_per_task $mem_per_task"
}
export -f get_slurm_specs

