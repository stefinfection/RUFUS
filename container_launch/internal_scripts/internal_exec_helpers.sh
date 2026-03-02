#!/bin/bash
# A collection of functions run inside the container to get command line arguments 
# for RUFUS and post-process execution.
# Requires the RUFUS config yaml to be parsed and translated into a cleaned env file.
# SJ Georges Nov2025

# ----------------- EXEC FUNCTIONS ----------------- #

# TODO: do I need to resource env + globals file here?

# Returns Fetches control or kg1 hash from S3 for region if region arg provided, or whole genome hash otherwise
# Returns path inside container to downloaded hash directory (named REMOTE_CONTROL/KG1_HASH_DIR in globals.env)
# Called indirectly from main runRufus.sh script (inside the container)
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
        bash "$DOWNLOAD_HASH_SCRIPT" "${hash_type}" "${hash_version}" "wg" "$remote_hash_dir" >&2
        hash="${remote_hash_dir}wg_${hash_type}_${hash_version}.Jhash"
    else
        # Convert chrN:n-m to chrN_n_m
        echo "Fetching version ${hash_version} ${hash_type} hash for region: $region" >&2
        local fmtd_reg=$(echo "$region" | tr ':-' '_')
        bash "${DOWNLOAD_HASH_SCRIPT}" "${hash_type}" "${hash_version}" "$fmtd_reg" "$remote_hash_dir" >&2
        hash="${remote_hash_dir}${fmtd_reg}_${hash_type}_${hash_version}.Jhash"
    fi

    echo "$hash"
}

# Looks for control or kg1 hashes (in $LOCAL_{CONTROL/KG1}_HASH_DIR in clean.env) first
# At this point, we know that if a local directory has been provided and mounted, it contains at least one *.Jhash file
# If local directory contains multiple *.Jhash files matching region or WG, will return error code
# If can't find, will pull from S3
# Called indirectly from main runRufus.sh script (inside container)
get_hash() {
    local region="$1"
    local geo_type="$2"
    local hash_type="$3"

    local hash_type_upper=$(echo "$hash_type" | tr '[:lower:]' '[:upper:]')

    # Output vars
    local hash=""
    local output_echo=""
    
    # Get the user's original directory path for error messages
    local external_local_hash_var="EXTERNAL_LOCAL_${hash_type_upper}_HASH_DIR" # Must correspond to that in rufus.env
    local external_local_hash_dir="${!external_local_hash_var}"

    # Get container visible directory path for actual work
    local internal_local_hash_var="INTERNAL_LOCAL_${hash_type_upper}_HASH_DIR" # Must correspond to that in globals.env
    local internal_local_hash_dir="${!internal_local_hash_var}"

    # Local hashes
    if [ "$geo_type" == "local" ]; then
        # Whole genome mode and we're looking locally
        if [ "$region" == "" ]; then
            # Check that we have a single wg file
            file_count=$(bash -c "ls ${internal_local_hash_dir}*wg*.Jhash 2>/dev/null | wc -l")    
            if [ "$file_count" -gt 1 ]; then
                echo "ERROR: Multiple whole genome ${hash_type} hash files found in ${external_local_hash_dir}:" >&2
                ls "${internal_local_hash_dir}"*wg*.Jhash >&2
                echo "Please ensure only one *wg*.Jhash file exists in the directory." >&2
                return 1
            elif [ "$file_count" -eq 1 ]; then
                hash=$(bash -c "ls ${internal_local_hash_dir}*wg*.Jhash")
                output_echo+="Using local whole genome ${hash_type}: $hash"
            else
                # File does not exist locally, fetch from S3
                output_echo+="Could not find local whole genome ${hash_type} hash in ${external_local_hash_dir}. The file in this directory must be named like \"*wg*.Jhash\" for RUFUS to recognize it. Fetching hash from RUFUS S3 repository."
                hash=$(fetch_hash "$region" "$hash_type")
            fi
        # Region mode and we're looking locally
        else
            fmtd_reg=$(echo "$region" | tr ':-' '_')
            file_count=$(bash -c "ls ${internal_local_hash_dir}*${fmtd_reg}*.Jhash 2>/dev/null | wc -l")
            if [ "$file_count" -gt 1 ]; then
                echo "ERROR: Multiple ${hash_type} hash files for region $region found in ${external_local_hash_dir}:" >&2
                ls "${internal_local_hash_dir}"*${fmtd_reg}*.Jhash >&2
                echo "Please ensure only one *${fmtd_reg}*.Jhash file exists in the directory." >&2
                return 1
            elif [ "$file_count" -eq 1 ]; then
                hash=$(bash -c "ls ${internal_local_hash_dir}*${fmtd_reg}*.Jhash")
                output_echo+="Using local ${hash_type} hash for region $region at: $hash"
            else
                # File does not exist locally, fetch from S3
                output_echo+="Could not find local ${hash_type} hash for region $region in ${external_local_hash_dir}. The file in this directory must be named like \"*${fmtd_reg}*.Jhash\" for RUFUS to recognize it. Fetching hash from RUFUS S3 repository."
                hash=$(fetch_hash "$region" "$hash_type")
            fi
        fi
    # Remote fetching of hashes
    else
        output_echo+="Fetching ${hash_type} hash from RUFUS S3 repository."
        hash=$(fetch_hash "$region" "$hash_type")
    fi
    
    echo -e "$output_echo" >&2
    echo "$hash"
}

# Returns RUFUS argument string for 1000G hashes as appropriate (without -e)
# Called by main rufus script
get_control_arg() {
    flag="$1"
    ctrl_arg=""

    if [ "$flag" == "local" ]; then
        ctrl_hash=$(get_hash $REGION "local" "control") || exit 1
    else
        ctrl_hash=$(get_hash $REGION "remote" "control") || exit 1
    fi
    ctrl_arg+="$ctrl_hash "

    echo "$ctrl_arg"
}

# Returns RUFUS argument string for 1000G hashes as appropriate (without -e)
# Called by main rufus script
get_kg1_arg() {
    flag="$1"
    kg1_arg=""

    if [ "$flag" == "local" ]; then
        kg1_hash=$(get_hash $REGION "local" "kg1") || exit 1
    else
        kg1_hash=$(get_hash $REGION "remote" "kg1") || exit 1
    fi
    kg1_arg+="$kg1_hash "

    echo "$kg1_arg"
}


# Returns RUFUS flag for control prebuilt hashes, if optioned
# Otherwise, returns empty string
# Run at setup time
get_control_hash_flag() {
    local ctrl_hash_arg=""

    # Get control region specific line, if optioned
    if [ "$NO_CONTROL_REMOVAL" != "TRUE" ]; then
        if [ -n "$EXTERNAL_LOCAL_CONTROL_HASH_DIR" ]; then
            ctrl_hash_arg+="--local-ctrl"
        else
            ctrl_hash_arg+="--remote-ctrl"
        fi
    fi

    echo "$ctrl_hash_arg"
}

# Returns RUFUS flag for kg1 prebuilt hashes, if optioned
# Otherwise, returns empty string
# Run at setup time
get_kg1_hash_flag() {
    local kg1_hash_arg=""

    # Get KG1 region specific line, if optioned
    if [ "$NO_KG1_REMOVAL" != "TRUE" ]; then
        if [ -n "$EXTERNAL_LOCAL_KG1_HASH_DIR" ]; then
            kg1_hash_arg+="--local-kg1"
        else
            kg1_hash_arg+="--remote-kg1"
        fi
    fi

    echo "$kg1_hash_arg"
}

# Returns RUFUS argument string for reference
# Agnostic to BWA index status
# Run at setup time
get_ref_arg() {
    ref_base=$(basename "$REFERENCE_FASTA")
    ref_arg="-r ${REF_INDEX_DIR}${ref_base}"

    # If subject_file ends with cram, need to change region_arg to -cr
    subject_base=$(basename "$SUBJECT_FILE")
    if [[ "$subject_base" == *.cram ]]; then
        ref_arg="-cr ${REF_INDEX_DIR}${ref_base}"
    fi
}

# Returns single string of arguments provided directly to RUFUS run script
# All args here are not relative to a region
# TODO: do I want to return array of strings here and do printf at PoC to ensure proper formatting?
get_rufus_args() {
    local rufus_args=""

    subject_base=$(basename "${SUBJECT_FILE}")

    if [ "${#CONTROL_FILE_ARRAY[@]}" -ne 0 ]; then
        # Concatenate controls into -c delimited string
        for control in "${CONTROL_FILE_ARRAY[@]}"; do
            control_base=$(basename "$control")
            ctrl_arg+="-c ${RUNTIME_TEMP_DIR}${control_base} "
        done
    fi

    ctrl_hash_flag=$(get_control_hash_flag)
    kg1_hash_flag=$(get_kg1_hash_flag)
    ref_arg=$(get_ref_arg)

    rufus_args+="\
    -s ${RUNTIME_TEMP_DIR}${subject_base} \
    $ctrl_arg \
    $ref_arg \
    $kg1_hash_flag \
    $ctrl_hash_flag \
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