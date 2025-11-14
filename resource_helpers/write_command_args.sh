#!/bin/bash
# This is run within container

# Constants
CMD_OUT="/mnt/rufus_supplementals/rufus.cmd"
GLOBALS_FILE="/opt/RUFUS/resources/globals.txt"

# Required args
CONTAINER_ID="$1"
CONT_ENV_FILE="/mnt/rufus.env"  # Mounted into container

# Make env variables available
if [ -f "$CONT_ENV_FILE" ]; then
    set -a
    source <(grep -v '^#' $CONT_ENV_FILE | grep -v '^[[:space:]]*$' | sed 's/\r$//')
    set +a
else
    echo "Error: $CONT_ENV_FILE file not found - please provide valid path to rufus.env file"
    exit 1
fi

# Make globals available
if [ -f "$GLOBALS_FILE" ]; then
    set -a
    source <(grep -v '^#' $GLOBALS_FILE | grep -v '^[[:space:]]*$' | sed 's/\r$//')
    set +a
else
    echo "Error: $GLOBALS_FILE file not found - please provide valid path to rufus.env file"
    exit 1
fi

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

run_cmd="docker exec ${CONTAINER_ID} bash /opt/RUFUS/runRufus.sh -s /mnt/$SUBJECT_FILE $ctrl_arg -r $REFERENCE_FASTA -k $KMER_LENGTH -m $KMER_DEPTH_CUTOFF -t $THREAD_LIMIT $OTHER_FLAGS -e kg1_$KG1_HASH_VERSION"
post_cmd="docker exec ${CONTAINER_ID} bash /opt/RUFUS/post_process/post_process.sh -s /mnt/$SUBJECT_FILE -r $REFERENCE_FASTA -w $WINDOW_SIZE -d /mnt $ctrl_arg"

mkdir -p /mnt/rufus_supplementals

echo "##RUFUSCommandLine=<ID=rufus, Branch=\"$RUFUS_BRANCH\", Version=\"$RUFUS_VERSION\", Command=\"$run_cmd\">" > "$CMD_OUT"
echo "##RUFUSCommandLine=<ID=rufus, Branch=\"$RUFUS_BRANCH\", Version=\"$RUFUS_VERSION\", Command=\"$post_cmd\">" >> "$CMD_OUT"
