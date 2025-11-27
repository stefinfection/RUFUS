#!/bin/bash
# A collection of functions run inside the container to get mount clauses for container launch.
# Requires the RUFUS config yaml to be parsed and translated into a cleaned env file.
# SJ Georges Nov2025

# TODO: do I need to resource here? everything in here called from build launch script?
# Import cleaned env file
# ENV_FILE="/mnt/rufus_temp/cleaned.env"
# set -a
# source <(grep -v '^#' $ENV_FILE | grep -v '^[[:space:]]*$' | sed 's/\r$//')
# set +a

# ----------------- MOUNT FUNCTIONS ----------------- #

# Returns array of subject file + index
get_subject_array() {
    local subj_array=()

    # Add index
    if [[ "$SUBJECT_FILE" == *.cram ]]; then
        post_fix="crai"
    elif [[ "$SUBJECT_FILE" == *.bam ]]; then
        post_fix="bai"
    fi

    subj_array=(
        "$SUBJECT_FILE"
        "$SUBJECT_FILE.$post_fix"
    )

    echo "${subj_array[@]}" #TODO: is this the right way to return array?
}
export -f get_subject_array

# Returns array of all paired controls + indexes, if they exist
get_paired_ctrl_array() {
    local ctrl_array=()

    if [ "${#CONTROL_FILE_ARRAY[@]}" -ne 0 ]; then
        # Concatenate controls into -c delimited string
        for control in "${CONTROL_FILE_ARRAY[@]}"; do
                if [[ "$control" == *.cram ]]; then
                    post_fix="crai"
                elif [[ "$control" == *.bam ]]; then
                    post_fix="bai"
                fi
                ctrl_array+=("$control" "$control.postfix")
        done
    fi

    echo "${ctrl_array[@]}" #TODO: is this the right way to return array?
}
export -f get_paired_ctrl_array

# Returns mount clause for local prebuilt hashes, if optioned
get_hash_mount() {
    local hash_type="$1"
    local mount_op="$2"

    local hash_clause=""

    if [ "$hash_type" == "ctrl" ] && [ -n "$EXTERNAL_LOCAL_CONTROL_HASH_DIR" ]; then
        control_path=$(dirname "$EXTERNAL_LOCAL_CONTROL_HASH_DIR")
        hash_clause="$mount_op $control_path:$INTERNAL_LOCAL_CONTROL_HASH_DIR"
    elif [ "$hash_type" == "kg1" ] && [ -n "$EXTERNAL_LOCAL_KG1_HASH_DIR" ]; then
        kg1_path=$(dirname "$EXTERNAL_LOCAL_KG1_HASH_DIR")
        hash_clause="$mount_op $kg1_path:$INTERNAL_LOCAL_KG1_HASH_DIR"
    fi

    echo "$hash_clause"
}
export -f get_hash_mount

# TODO: should I just mount to dir here instead of individually? Enforcing they have to be in same dir anyways
# Returns array of reference file + all corresponding BWA indexes
get_ref_array() {
    local ref_array=("$REFERENCE_FASTA")

    # Determine the base filename (without .gz if present)
    if [[ "$REFERENCE_FASTA" == *.gz ]]; then
        ref_file="${REFERENCE_FASTA%.gz}"
    else
        ref_file="$REFERENCE_FASTA"
    fi

    # Add indexes
    for ext in sa bwt pac amb ann fai; do
        ref_array+=("${ref_file}.${ext}")
    done

    echo "${ref_array[@]}" # TODO: is this the correct way to return array?
}
export -f get_ref_array

# Returns mount clause for controls, 1000G, subject, and references/indexes
# For either container technology
get_container_mount_clause() {
    local container_type="$1"

    # Get option flag based on container type
    mount_opt="-v"
    if [ "$container_type" == "$SINGULARITY" ]; then
        mount_opt="--bind"
    fi

    mount_lines=()

    # Add subject file + index
    subject_mount_array=$(get_subject_array)
    for sub in "${subject_mount_array[@]}"; do
        sub_basename=$(basename $sub)
        mount_lines+=("$mount_opt $sub:${RUNTIME_TEMP_DIR}${sub_basename}:ro")
    done

    # Add paired controls, if optioned
    control_mount_array=$(get_paired_ctrl_array)
    for ctrl in "${control_mount_array[@]}"; do
        ctrl_basename=$(basename "$ctrl")
        mount_lines+=("$mount_opt $ctrl:${RUNTIME_TEMP_DIR}${ctrl_basename}:ro")
    done

    # Add control hashes, if optioned
    ctrl_hash_clause=$(get_hash_mount "ctrl" "$mount_op")
    mount_lines+=("$ctrl_hash_clause")
    
    # Add 1000G hashes, if optioned
    kg1_hash_clause=$(get_hash_mount "kg1" "$mount_op")
    mount_lines+=("$kg1_hash_clause")

    # Add reference + bwa indexes
    ref_mount_array=$(get_ref_array)
    for ref in "${ref_mount_array[@]}"; do
        mount_lines+=("$mount_opt $ref:${REF_INDEX_DIR}")
    done

    # Write out readably
    mount_clause=""
    for clause in "${mount_lines[@]}"; do
        printf -v mount_clause '%s%s \\\n' "$mount_clause" "$clause"
    done

    echo "$mount_clause"
}
export -f get_container_mount_clause