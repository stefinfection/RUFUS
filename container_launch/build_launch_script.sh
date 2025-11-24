#!/bin/bash
# Generates and writes launch script according to args in rufus_config.yaml
# SJ Georges Nov2025

container_type="$1"
host_config_file="$2"

# Constants
config_file=/temp/rufus_config.yaml
SINGULARITY="singularity"
DOCKER="docker"

# Import helpers
source "/opt/RUFUS/container_launch/internal_launch_helpers.sh"


#-------------------WRITE FUNCTIONS -------------------#

# Writes out #SBATCH header or bash directive according to container type
# Writes out small instruction piece, date stamp, and version
write_header() {
    local out_script="$1"
    # get #SBATCH header file for job run script (TODO: how to code this human readable that will echo correctly

    if [ "$container_type" == "$SINGULARITY" ]; then

        # Output required slurm args
        cat <<EOF > "$out_script"
#!/bin/bash
#SBATCH --nodes=${SLURM_NODES}
#SBATCH --account=${SLURM_ACCOUNT}
#SBATCH --partition=${SLURM_PARTITION}
#SBATCH --time=${SLURM_TIME}
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err
EOF

        # Output optional slurm args
        if [ -n "${SLURM_JOB_NAME}" ]; then
            echo "#SBATCH --job-name=${SLURM_JOB_NAME}" >> "$out_script"
        else
            echo "#SBATCH --job-name=rufus" >> "$out_script"
        fi

        if [ -n "${SLURM_EMAIL}" ]; then
            echo "#SBATCH --mail-type=ALL" >> "$out_script"
            echo "#SBATCH --mail-user=${SLURM_EMAIL}" >> "$out_script"
        fi

        # Directions
        cat <<EOF >> "$out_script"

This script was automatically generated with RUFUS ${RUFUS_VERSION} based on the config file ${host_config_file}.
Review, and then launch with: sbatch "$out_script"

EOF

        # This function outputs the script that will be launched once on each node
    else
        echo -e "#!/bin/bash" > "$out_script"
        echo -e "This script was automatically generated with RUFUS ${RUFUS_VERSION} based on the config file ${host_config_file}." >> "$out_script"
        echo -e "Review, and then launch with: sbatch $out_script" >> "$out_script"
        echo "" >> "$out_script" # insert empty line
    fi
}
export -f write_header

# Writes out the script portion that binds volumes and starts container
# Writes out fail trap
write_container_start_piece() {
    local out_script="$1"

    if [ "$container_type" == "$SINGULARITY" ]; then

    else
        echo "$CONTAINER_ID=$(docker run -d --rm --name rufus-worker \
        -u "$(id -u):$(id -g)" \
        -v "${WORKING_DIR}:/mnt" \
        --cap-add SYS_ADMIN \ 
        --device /dev/fuse \" >> $out_script

        local mount_clause=$(get_container_mount_clause)
        echo "$mount_clause"
        
        echo "$RUFUS_DOCKER_IMAGE" \
            tail -f /dev/null)" >> $out_script
    fi

}

write_container_stop_piece() {
    local out_script="$1"

    if [ "$container_type" == "$SINGULARITY" ]; then

    else

    fi
}

write_bwa_indexes() {
    # TODO: wrap this in appropriate srun if singularity + depend regional jobs on finish
    echo "# Generate BWA indexes for reference fasta" >> "$launch_out"
    echo "# Note: this step can be skipped next time you run by copying indexes into same location of $REFERENCE_FASTA" >> "$launch_out"
    echo "docker exec ${CONTAINER_ID} bash /opt/RUFUS/resource_helpers/build_bwa_indexes.sh ${REFERENCE_FASTA}" >> "$launch_out"
}

write_rufus_execution_piece() {
    local out_script="$1"

    # pull out region worker/internal launch helpers script
    # start container with srun

    # If attempt to use more nodes than region (in case doing whole genome mode), warn only using one


    # write out regional specific arg fxn calls
    # It gets the regional specific arguments for control and kg1 hashes (these functions will be internal to rufus container now)
    # It sruns the individual region jobs
    # use srun to start array jobs, instead of writing another script to sbatch (just fill in cpus-per-task and mem to have slurm max parallelism)
    # try 8G mem and 10 cpus per task (so 20 threads to rufus) if region - otherwise do whole genome

}

write_post_process_piece() {
    local out_script="$1"

    if [ "$container_type" == "$SINGULARITY" ]; then

    else

    fi

    # write post-process script
    # slurm_id=sbatch slurm script
    # sbatch post-process script --depend-on:slurm_id
    # If first script in all nodes queued (FIRST_NODE=$(scontrol show hostnames "$SLURM_JOB_NODELIST" | head -n1)), then run post-process
    # get control argument formatted for post process (will use for both so do before if statement)

}

write_out_script() {
    if [ "$container_type" == "$SINGULARITY" ]; then
        out_script="${WORKING_DIR}/run_rufus_singularity.slurm"
    else
        out_script="${WORKING_DIR}/run_rufus_docker.sh"
    fi
        write_header

        write_container_start_piece

        echo "start_time=$(date +%s)" >> $out_script

        # write out arguments to rufus.cmd for final vcf header

        write_rufus_execution_piece

        # TODO: echo time stamp for rufus small jobs piece

        write_post_process_piece

        # TODO: echo total run time including post-processing

        write_container_stop_piece

}


#------------------- DO WORK -------------------#

# Parse config file and write to temp.env file
parse_config_file

# Write BWA 
exist=$(check_for_bwa_indexes)
if [ ${exist} == "false" ]; then
    write_bwa_indexes
fi

write_out_script
