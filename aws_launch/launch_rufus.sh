#!/bin/bash

# CONSTANTS
#DEV_MOUNT="-v /home/ubuntu/RUFUS/runRufus.sh:/opt/RUFUS/runRufus.sh \
#  -v /home/ubuntu/RUFUS/scripts:/opt/RUFUS/scripts \
#  -v /home/ubuntu/RUFUS/resource_helpers:/opt/RUFUS/resource_helpers \
#  -v /home/ubuntu/RUFUS/post_process:/opt/RUFUS/post_process \
#  -v /home/ubuntu/RUFUS/resources:/opt/RUFUS/resources
#  -v /home/ubuntu/RUFUS/aws_launch/process_region_worker.sh:/opt/RUFUS/aws_launch/process_region_worker.sh \
#  -v /home/ubuntu/RUFUS/bin/RUFUS.interpret:/opt/RUFUS/bin/RUFUS.interpret"
DEV_MOUNT=""

# Check for required argument
ENV_FILE="$1"
if [ -z "$ENV_FILE" ]; then
    echo "Error: Please provide PATH_TO_RUFUS_ENV argument" >&2
    exit 1
fi

# Check RUFUS env file arg actually exists
if [ -f "$ENV_FILE" ]; then
    ENV_FILE=$(realpath "$ENV_FILE")
    set -a
    source <(grep -v '^#' $ENV_FILE | grep -v '^[[:space:]]*$' | sed 's/\r$//')
    set +a
else
    echo "Error: $ENV_FILE file not found - please provide valid path to rufus.env file"
    exit 1
fi

# If we don't have a working dir, set it to .
if [ -z "$WORKING_DIR" ]; then
    WORKING_DIR=$(pwd)
    WORKING_DIR=$(realpath "$WORKING_DIR")
fi

# Make temp env file with realpaths for all input files
TEMP_ENV_FILE=${WORKING_DIR}/temp_rufus.env
touch $TEMP_ENV_FILE
cat $ENV_FILE > $TEMP_ENV_FILE

# Checks for required arguments and formatting + bounds of integer arguments
check_inputs() {
    # Check for all required variables to be filled in rufus.env
    REQUIRED_VARS=(SUBJECT_FILE KMER_DEPTH_CUTOFF THREAD_LIMIT REFERENCE_FASTA RUFUS_DOCKER_IMAGE)
    for var in "${REQUIRED_VARS[@]}"; do
        if [ -z "${!var}" ]; then
            echo "Error: Required variable $var is not set in rufus.env"
            exit 1
        fi
    done

    # Check for subject file existence
    subject_path=$(realpath "${SUBJECT_FILE}")
    if [ ! -f "$subject_path" ]; then
        echo "Error: SUBJECT_FILE $subject_path not found" >&2
        exit 1
    else
        echo "SUBJECT_FILE=$subject_path" >> $TEMP_ENV_FILE
    fi

    # Check for control files existence
    temp_ctrl_array=()
    for control in "${CONTROL_FILE_ARRAY[@]}"; do
        control_path=$(realpath "$control")
        if [ ! -f "$control_path" ]; then
            echo "Error: Control file $control_path not found" >&2
            exit 1
        else
            temp_ctrl_array+=("$control_path")
        fi
    done
    CONTROL_FILE_ARRAY=("${temp_ctrl_array[@]}")
    echo "CONTROL_FILE_ARRAY=($temp_ctrl_array)" >> $TEMP_ENV_FILE

    # Check for reference fasta existence
    reference_path=$(realpath "${REFERENCE_FASTA}")
    if [ ! -f "$reference_path" ]; then
        echo "Error: REFERENCE_FASTA $reference_path not found" >&2
        exit 1
    else
	REFERENCE_FASTA=$reference_path
        echo "REFERENCE_FASTA=$reference_path" >> $TEMP_ENV_FILE
    fi

    # Check if region file provided, that it exists
    if [ -n "$REGION_FILE" ]; then
        region_path=$(realpath "${REGION_FILE}")
        if [ ! -f "$region_path" ]; then
            echo "Error: REGION_FILE $region_path not found. Please provide valid file or leave empty for whole genome mode." >&2
            exit 1
        else
	    REGION_FILE=$region_path
            echo "REGION_FILE=$region_path" >> $TEMP_ENV_FILE
        fi
    fi

    # Check that thread limit is an int greater than 0
    if ! [[ "$THREAD_LIMIT" =~ ^[0-9]+$ ]] || [ "$THREAD_LIMIT" -le 0 ]; then
        echo "Error: THREAD_LIMIT must be a positive integer" >&2
        exit 1
    fi


    # Check that if one of the following are filled out, the other two also are - WINDOW_SIZE, JOB_THRESHOLD, REGION_FILE
    if { [ -n "$WINDOW_SIZE" ] || [ -n "$JOB_THRESHOLD" ] || [ -n "$REGION_FILE" ]; } && { [ -z "$WINDOW_SIZE" ] || [ -z "$JOB_THRESHOLD" ] || [ -z "$REGION_FILE" ]; }; then
        echo "Error: If one of WINDOW_SIZE, JOB_THRESHOLD, or REGION_FILE is filled out, all three must be provided for regional processing mode" >&2
        exit 1
    fi

    # Check that if job threshold is present, it is an int greater than 0
    if [ -n "$JOB_THRESHOLD" ]; then
        if ! [[ "$JOB_THRESHOLD" =~ ^[0-9]+$ ]] || [ "$JOB_THRESHOLD" -le 0 ]; then
            echo "Error: JOB_THRESHOLD must be a positive integer" >&2
            exit 1
        fi
    fi

    # Check that kmer depth cutoff is an int and warn if less than 3
    if ! [[ "$KMER_DEPTH_CUTOFF" =~ ^[0-9]+$ ]]; then
        echo "Error: KMER_DEPTH_CUTOFF must be a positive integer" >&2
        exit 1
    elif [ "$KMER_DEPTH_CUTOFF" -lt 3 ]; then
        echo "Warning: KMER_DEPTH_CUTOFF is set to less than 3, which may lead to increased false positives" >&2
    fi

    # Check if KMER_LENGTH filled out, is an int and warn if not 25
    if [ -n "$KMER_LENGTH" ]; then
        if ! [[ "$KMER_LENGTH" =~ ^[0-9]+$ ]]; then
            echo "Error: KMER_LENGTH must be a positive integer" >&2
            exit 1
        elif [ "$KMER_LENGTH" -ne 25 ]; then
            echo "Warning: KMER_LENGTH is set to $KMER_LENGTH, RUFUS has been robustly tested with a KMER_LENGTH of 25 and is recommended" >&2
        fi
    else
        KMER_LENGTH=25
        echo "KMER_LENGTH=$KMER_LENGTH" >> $TEMP_ENV_FILE 
    fi

    # Check if WINDOW_SIZE filled out, if it is, make sure an int and is 1000
    if [ -n "$WINDOW_SIZE" ]; then
        if ! [[ "$WINDOW_SIZE" =~ ^[0-9]+$ ]]; then
            echo "Error: WINDOW_SIZE must be a positive integer" >&2
            exit 1
        elif [ "$WINDOW_SIZE" -ne 1000 ]; then
            echo "ERROR: WINDOW_SIZE is set to $WINDOW_SIZE - must be either 1000 or empty" >&2
            exit 1
        fi
    fi

    echo "WORKING_DIR=$WORKING_DIR" >> $TEMP_ENV_FILE
}
export -f check_inputs

# Check for correct controls setup and returns paths needed for mounting if necessary
# Sets up link to realpath of file within provided directory, because may be a symlink
# WARNING: all Jhash files must be in the same realpath directory for mounting to work correctly
# Hashes get mounted to /mnt/rufus_resources/control_hashes within container
# Paired controls get mounted to /mnt/paired_controls/
set_up_controls() {
    # Check for controls here and notify if using internal
    if [ ${#CONTROL_FILE_ARRAY[@]} -eq 0 ]; then
        echo "No paired controls provided, running RUFUS in internal control mode..." >&2
    fi

    local mount_clause=""
    # Make sure the local directory exists if provided
    if [ ! -z "${CONTROL_HASH_LOCAL_DIR}" ]; then
        control_path=$(realpath "${CONTROL_HASH_LOCAL_DIR}")
        if [ ! -d "${control_path}" ]; then
            echo "Error: CONTROL_HASH_LOCAL_DIR ${control_path} does not exist. Please ensure directory exists or leave CONTROL_HASH_LOCAL_DIR empty for S3 fetching." >&2
            return 1
        else
            # Make sure directory has at least one *.Jhash file in it
            shopt -s nullglob
            jhash_files=("${control_path}"/*.Jhash)
            shopt -u nullglob

            if [ ${#jhash_files[@]} -eq 0 ]; then
                echo "Error: CONTROL_HASH_LOCAL_DIR ${control_path} does not contain any *.Jhash files. Please ensure directory has Jhash files or leave CONTROL_HASH_LOCAL_DIR empty for S3 fetching." >&2
                return 1
            fi

            # Verify at least one symlink target is accessible
            accessible=false
            for file in "${jhash_files[@]}"; do
                if [ -f "$file" ]; then
                    file_path=$(realpath "$file")
                    parent_dir_file=$(dirname "$file_path")
                    # Mount to realpath of file rather than parent dir because file may be symlinked
                    mount_clause="-v ${parent_dir_file}:/mnt/rufus_resources/control_hashes "
                    echo "CONTROL_HASH_LOCAL_DIR=${parent_dir_file}" >> $TEMP_ENV_FILE
                    accessible=true
                    break
                fi
            done

            if [ "$accessible" = false ]; then
                echo "Error: *.Jhash files found but none are accessible (broken symlinks or goofys issue)." >&2
                return 1
            fi
        fi
    fi

    # If we have paired controls provided, also use those
    if [ "${#CONTROL_FILE_ARRAY[@]}" -ne 0 ]; then
        # Concatenate controls into -c delimited string
        ctrl_arg=""
        for control in "${CONTROL_FILE_ARRAY[@]}"; do
            ctrl_path=$(realpath "$control")
            ctrl_basename=$(basename "$ctrl_path")
            ctrl_arg+="-v $ctrl_path:/mnt/paired_controls/$ctrl_basename "
            
            # Add index or error if not found
            if [[ "$ctrl_basename" == *.bam ]]; then
                if [ -f "${ctrl_path}.bai" ]; then
                    ctrl_arg+="-v ${ctrl_path}.bai:/mnt/paired_controls/${ctrl_basename}.bai "
                elif [ -f "${ctrl_path%.bam}.bai" ]; then
                    bai_path="${ctrl_path%.bam}.bai"
                    ctrl_arg+="-v $bai_path:/mnt/paired_controls/$(basename $bai_path) "
                else
                    echo "ERROR: Could not find index file for control BAM ${ctrl_path}. Please ensure .bai file exists." >&2
                    return 1
                fi
            elif [[ "$ctrl_basename" == *.cram ]]; then
                if [ -f "${ctrl_path}.crai" ]; then
                    ctrl_arg+="-v ${ctrl_path}.crai:/mnt/paired_controls/${ctrl_basename}.crai "
                elif [ -f "${ctrl_path%.cram}.crai" ]; then
                    crai_path="${ctrl_path%.cram}.crai"
                    ctrl_arg+="-v $crai_path:/mnt/paired_controls/$(basename $crai_path) "
                else
                    echo "ERROR: Could not find index file for control CRAM ${ctrl_path}. Please ensure .crai file exists." >&2
                    return 1
                fi
            fi
        done
        mount_clause+="$ctrl_arg"
    fi
    echo "$mount_clause"
}
export -f set_up_controls

# Check for correct 1000G setup and returns paths needed for mounting if necessary
# Sets up link to realpath of file within provided directory, because may be a symlink
# WARNING: all Jhash files must be in the same realpath directory for mounting to work correctly
# Hashes get mounted to /mnt/rufus_resources/kg1_hashes within container
set_up_kg1() {
    local mount_clause=""

    # Notify if not removing 1000G variants
    if [ "$NO_KG1_REMOVAL" == "true" ] || [ "$NO_KG1_REMOVAL" == "TRUE" ]; then
        echo "Warning: not removing common population variants in the 1000G cohort" >&2
    else
        # Make sure the local directory exists if provided
        if [ ! -z "${KG1_HASH_LOCAL_DIR}" ]; then
            kg1_path=$(realpath "${KG1_HASH_LOCAL_DIR}")
            if [ ! -d "${kg1_path}" ]; then
                echo "Error: KG1_HASH_LOCAL_DIR ${kg1_path} does not exist. Please ensure directory exists or leave KG1_HASH_LOCAL_DIR empty for S3 fetching." >&2
                return 1
            else
                # Make sure directory has at least one *.Jhash file in it
                shopt -s nullglob
                jhash_files=("${kg1_path}"/*.Jhash)
                shopt -u nullglob
                if [ ${#jhash_files[@]} -eq 0 ]; then
                    echo "Error: KG1_HASH_LOCAL_DIR ${kg1_path} does not contain any *.Jhash files. Please ensure directory has Jhash files or leave KG1_HASH_LOCAL_DIR empty for S3 fetching." >&2
                    return 1
                fi

                # Verify at least one symlink target is accessible
                accessible=false
                for file in "${jhash_files[@]}"; do
                    if [ -f "$file" ]; then
                        file_path=$(realpath "$file")
                        parent_dir_file=$(dirname "$file_path")
                        # Mount to realpath of file rather than parent dir because file may be symlinked
                        mount_clause="-v ${parent_dir_file}:/mnt/rufus_resources/kg1_hashes "
                        echo "KG1_HASH_LOCAL_DIR=${parent_dir_file}" >> $TEMP_ENV_FILE
                        accessible=true
                        break
                    fi
                done

                if [ "$accessible" = false ]; then
                    echo "Error: *.Jhash files found but none are accessible (broken symlinks or goofys issue)." >&2
                    return 1
                fi
            fi
        fi
    fi
    echo "$mount_clause"
}
export -f set_up_kg1

set_up_ref() {

    local ref_path=$(realpath "${REFERENCE_FASTA}") # Absolute path on host machine
    local ref_base=$(basename ${ref_path})          # Base filename
    local path_to_ref="$(dirname ${ref_path})"      # Directory on host machine

    local build_refs="FALSE"

    # We'll assume indexes are in same dir as reference unless found otherwise
    local mount_clause="-v ${path_to_ref}:/mnt/bwa_indexes "
        
    # Determine the base filename (without .gz if present)
    if [[ "$ref_path" == *.gz ]]; then
        ref_file="${ref_path%.gz}"
        ref_file_base="$(basename ${ref_file})"
    else
        ref_file="$ref_path"
        ref_file_base="$ref_base"
    fi
    
    # Check all required index files in one loop
    for ext in sa bwt pac amb ann fai; do
	echo "checking for ${ref_file}.${ext}" >&2
        if [[ ! -e "${ref_file}.${ext}" ]]; then
            build_refs="TRUE"
            break
        fi
    done
    
    # Check if they exist in bwa_indexes sub-directory if not found above
    if [ "$build_refs" == "TRUE" ]; then
        build_refs="FALSE"
        for ext in sa bwt pac amb ann fai; do
            if [[ ! -e "${path_to_ref}/bwa_indexes/${ref_file_base}.${ext}" ]]; then
                build_refs="TRUE"
                break
            fi
        done
        
        if [ "$build_refs" == "FALSE" ]; then
            mount_clause="-v ${path_to_ref}/bwa_indexes:/mnt/bwa_indexes "
        fi
    fi
    
    echo "$mount_clause|$build_refs"
}
export -f set_up_ref

# Start work
check_inputs
IFS='|' read -r ref_mount build_refs < <(set_up_ref)
control_mount=$(set_up_controls) || exit 1
kg1_mount=$(set_up_kg1) || exit 1
subject_path=$(realpath "${SUBJECT_FILE}")
subject_base=$(basename ${subject_path})
input_mount_clause="$ref_mount $control_mount $kg1_mount -v ${subject_path}:/mnt/${subject_base} "

# Check for subject index if bam or cram
if [[ "$subject_base" == *.bam ]]; then
    if [ -f "${subject_path}.bai" ]; then
        input_mount_clause+="-v ${subject_path}.bai:/mnt/${subject_base}.bai "
    elif [ -f "${subject_path%.bam}.bai" ]; then
        bai_path="${subject_path%.bam}.bai"
        input_mount_clause+="-v $bai_path:/mnt/$(basename $bai_path) "
    else
        echo "ERROR: Could not find index file for subject BAM ${subject_path}. Please ensure .bai file exists." >&2
        exit 1
    fi
elif [[ "$subject_base" == *.cram ]]; then
    if [ -f "${subject_path}.crai" ]; then
        input_mount_clause+="-v ${subject_path}.crai:/mnt/${subject_base}.crai "
    elif [ -f "${subject_path%.cram}.crai" ]; then
        crai_path="${subject_path%.cram}.crai"
        input_mount_clause+="-v $crai_path:/mnt/$(basename $crai_path) "
    else
        echo "ERROR: Could not find index file for subject CRAM ${subject_path}. Please ensure .crai file exists." >&2
        exit 1
    fi
fi

CONTAINER_ID=$(docker run -d --rm --name rufus-worker \
  -v ${WORKING_DIR}:/mnt \
  --cap-add SYS_ADMIN \
  --device /dev/fuse \
  $input_mount_clause \
  $DEV_MOUNT \
  $RUFUS_DOCKER_IMAGE \
  tail -f /dev/null)

start_time=$(date +%s)

# Write commands for final vcf (do inside container so have access to RUFUS versioning)
docker exec ${CONTAINER_ID} bash /opt/RUFUS/resource_helpers/write_command_args.sh "$CONTAINER_ID" "$TEMP_ENV_FILE"

# Pull out worker script
#docker cp ${CONTAINER_ID}:/opt/RUFUS/aws_launch/process_region_worker.sh ${WORKING_DIR}/process_region_worker.sh
PR_WORKER="${WORKING_DIR}/process_region_worker_temp.sh"

# Check for BWA indexes and create if necessary
if [ "$build_refs" == "TRUE" ]; then
    echo "Generating BWA indexes for reference fasta..."
    docker exec ${CONTAINER_ID} bash /opt/RUFUS/resource_helpers/build_bwa_indexes.sh "${REFERENCE_FASTA}"
fi

# Make rufus resource dirs inside container
docker exec ${CONTAINER_ID} mkdir -p /mnt/rufus_resources/control_hashes
docker exec ${CONTAINER_ID} mkdir -p /mnt/rufus_resources/kg1_hashes
docker exec ${CONTAINER_ID} mkdir -p /mnt/rufus_resources/logs

# Start work
echo "Starting RUFUS job(s)..."
if [ -n "$REGION_FILE" ]; then
    parallel -j "$JOB_THRESHOLD" bash ${PR_WORKER} "$CONTAINER_ID" "$TEMP_ENV_FILE" {} :::: "$REGION_FILE"
else
    bash ${WORKING_DIR}/process_region_worker.sh "$CONTAINER_ID" "$TEMP_ENV_FILE" ""
fi

concat_ctrl_post_arg=""
if [ ${#CONTROL_FILE_ARRAY[@]} -gt 0 ]; then
    # Concatenate controls without -c delimiters
    concat_ctrls=""
    for control in "${CONTROL_FILE_ARRAY[@]}"; do
        ctrl_base=$(basename ${control})
        concat_ctrls+="/mnt/$ctrl_base "
    done
    concat_ctrl_post_arg="-c $concat_ctrls"
fi

subject_base=$(basename ${SUBJECT_FILE})
ref_base=$(basename ${REFERENCE_FASTA})

# Wait for all jobs to finish before combining + post-processing
echo "All RUFUS regional jobs completed. Starting merge and post-process..."
docker exec ${CONTAINER_ID} bash /opt/RUFUS/post_process/post_process.sh -s "/mnt/$subject_base" -r "/mnt/${ref_base}" -w "$WINDOW_SIZE" -d "/mnt" "$concat_ctrl_post_arg"

# Move any files we want to keep into working dir
keep_items=("/mnt/rufus_resources/logs")
for item in $keep_items; do
    mv $item $WORKING_DIR
done

# Stop container and clean up
docker stop rufus-worker
rm ${WORKING_DIR}/process_region_worker.sh
rm $TEMP_ENV_FILE

end_time=$(date +%s)
elapsed=$((end_time - start_time))
echo "RUFUS completed. Total run time: $elapsed"

