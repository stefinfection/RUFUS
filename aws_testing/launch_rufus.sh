#!/bin/bash
DEV_MOUNT="-v /home/ubuntu/RUFUS/runRufus.sh:/opt/RUFUS/runRufus.sh \
  -v /home/ubuntu/RUFUS/scripts:/opt/RUFUS/scripts \
  -v /home/ubuntu/RUFUS/resource_helpers:/opt/RUFUS/resource_helpers \
  -v /home/ubuntu/RUFUS/post_process:/opt/RUFUS/post_process \
  -v /home/ubuntu/RUFUS/resources:/opt/RUFUS/resources
  -v /home/ubuntu/RUFUS/bin/RUFUS.interpret:/opt/RUFUS/bin/RUFUS.interpret"

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

# Check for regions file to determine if running in WG mode or region mode
REGION_PATH=""
if [ -z "$REGION_FILE" ]; then
    REGION_PATH="${HOST_DATA_DIR}/rufus_resources/empty_region.txt"
    touch $REGION_PATH
    echo "Running RUFUS in whole genome mode"
else
    REGION_PATH="${HOST_DATA_DIR}/rufus_resources/$REGION_FILE"
    if [ -f "$REGION_PATH" ]; then
        num=$(cat "$REGION_PATH" | wc -l)
        first=$(cat "$REGION_PATH" | head -n 1)
        last=$(cat "$REGION_PATH" | tail -n 1)
        echo "Using $REGION_PATH with $num regions - first is $first and last is $last"
    else
        echo "Error: REGION_FILE $REGION_PATH not found" >&2
        exit 1
    fi
fi

# Check for process_region_worker.sh script
process_region_worker="${HOST_DATA_DIR}/rufus_resources/process_region_worker.sh"
if [ ! -f "$process_region_worker" ]; then
    echo "Error: /process_region_worker.sh script not found in container" >&2
    exit 1
else
    chmod +x "$process_region_worker"
fi

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

# Returns docker command line invocation for putting in final vcf
get_invocation() {
    type="$1"
    container_id="$2"

    # get control argument or internal version
    ctrl_arg=""
    if [ "${#CONTROL_FILE_ARRAY[@]}" -eq 0 ]; then
        ctrl_arg="-e internal_$CONTROL_HASH_VERSION"
    else
        # Concatenate controls into a single -c delimited string
        for control in "${CONTROL_FILE_ARRAY[@]}"; do
            ctrl_arg+="-c /mnt/$control"
        done
    fi

    if [ "$type" == "post_process" ]; then
        echo "docker exec ${container_id} bash /opt/RUFUS/post_process/post_process.sh -s \"/mnt/$SUBJECT_FILE\" -r \"$(get_reference)\" -w \"$WINDOW_SIZE\" -d \"/mnt\" $ctrl_arg"
    else
        echo "docker exec ${container_id} bash /opt/RUFUS/runRufus.sh -s \"/mnt/$SUBJECT_FILE\" $ctrl_arg -r \"$(get_reference)\" -k \"$KMER_LENGTH\" -m \"$KMER_DEPTH_CUTOFF\" -t \"$THREAD_LIMIT\" $OTHER_FLAGS -e kg1_$KG1_HASH_VERSION"
    fi
}
export -f get_invocation

# Start RUFUS container
CONTAINER_ID=$(docker run -d --rm --name rufus-worker \
  -v /mnt/data:/mnt \
  $DEV_MOUNT \
  rufus:latest \
  tail -f /dev/null)

# Check for controls here and notify if using internal
if [ ${#CONTROLS[@]} -eq 0 ]; then
    echo "No controls provided, running RUFUS in internal control mode..."
fi

# TODO: left off here - test this
# Write commands for final vcf
cmd="$(get_invocation run_rufus $CONTAINER_ID)"
echo -e $cmd > "${HOST_DATA_DIR}/rufus_resources/rufus.cmd"

cmd="$(get_invocation post_process $CONTAINER_ID)"
echo -e $cmd > "${HOST_DATA_DIR}/rufus_resources/rufus.cmd"

# Start work
echo "Starting RUFUS job(s)..."
start_time=$(date +%s)

parallel -j "$JOB_THRESHOLD" "${process_region_worker}" "$ENV_FILE" "$CONTAINER_ID" {} :::: "$REGION_PATH"

concat_ctrl_post_arg=""
if [ ${#CONTROLS[@]} -gt 0 ]; then
    # Concatenate controls without -c delimiters
    concat_ctrls=""
    for control in "${CONTROL_FILE_ARRAY[@]}"; do
        concat_ctrls+="/mnt/$control "
    done
    concat_ctrl_post_arg="-c $concat_ctrls"
fi

ref=$(get_reference)

# Wait for all jobs to finish before combining + post-processing
echo "All RUFUS regional jobs completed. Starting merge and post-process..."
docker exec ${CONTAINER_ID} bash /opt/RUFUS/post_process/post_process.sh -s "/mnt/$SUBJECT_FILE" -r "$ref" -w "$WINDOW_SIZE" -d "/mnt" "$concat_ctrl_post_arg"

# Stop container and clean up
docker stop rufus-worker

end_time=$(date +%s)
elapsed=$((end_time - start_time))
echo "RUFUS completed. Total run time: $elapsed"
