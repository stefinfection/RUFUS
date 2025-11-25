#!/bin/bash
# Generates and writes launch script according to args in rufus_config.yaml
# Run from within the container.
# SJ Georges Nov2025

container_type="$1"
host_config_file="$2"

# Import cleaned env file here to access globals
ENV_FILE="/mnt/rufus_temp/cleaned.env"
set -a
source <(grep -v '^#' $ENV_FILE | grep -v '^[[:space:]]*$' | sed 's/\r$//')
set +a

# Import helper functions
source "${MOUNT_HELPERS}"
source "${EXEC_HELPERS}"

# Constants
SINGULARITY="singularity"
PR_WRAPPER="${WORKING_DIR}rufus_temp/process_region.sh" #TODO: need to pull this out and put in a temp dir for host access if using parallel

# ------------------- WRITE FUNCTIONS ------------------- #

# Writes out #SBATCH header or bash directive according to container type
# Writes out small instruction piece and version
write_header() {
    local out_script="$1"

    # Singularity version run on SLURM_NODES nodes
    if [ "$container_type" == "$SINGULARITY" ]; then

        # Output required slurm args
        cat <<EOF > "$out_script"
#!/bin/bash
#SBATCH --nodes=${SLURM_NODES}
#SBATCH --account=${SLURM_ACCOUNT}
#SBATCH --partition=${SLURM_PARTITION}
#SBATCH --time=${SLURM_TIME}
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err"
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
Review, and then launch with: sbatch $out_script

EOF


    # Docker version run on a single node
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
    local mount_clause=$(get_container_mount_clause "$container_type")
    echo "Start single container per node..." >> "$out_script"

    
    if [ "$container_type" == "$SINGULARITY" ]; then
        
        # Write out singularity trap
        cat <<EOF >> "$out_script"
cleanup() {
    echo "Something went wrong, stopping instances on all nodes..."
    # attempt to stop on every allocated node
    srun --ntasks-per-node=1 --exclusive \
        bash -lc "singularity instance stop rufus_instance_\${SLURM_NODEID} || true"
}
trap cleanup EXIT

EOF

        # Start one instance per node (blocking: returns after instance started on each node)
        cat <<EOF >> "$out_script"
# Start one instance per node
srun --ntasks-per-node --exclusive
bash -lc \"hostname; singularity instance start \
--writable-tmpfs \
--bind "${WORKING_DIR}:/work" \
"$mount_clause"
"$RUFUS_SINGULARITY_IMAGE" \
"rufus_instance_\${SLURM_NODEID}"\"

EOF

    else

        # Write out docker trap
        cat <<EOF >> "$out_script"
cleanup() {
    echo "Something went wrong, stopping container..."
    if docker ps -a --format '{{.Names}}' | grep -q "^rufus_worker\$"; then
        docker stop rufus_worker >/dev/null 2>&1 || true
    fi
}
trap cleanup EXIT

EOF
        # Start docker container on single node
        cat <<EOF >> "$out_script" 
\${CONTAINER_ID}=\$(docker run -d --rm --name rufus-worker \
-u "\$(id -u):\$(id -g)" \
-v "${WORKING_DIR}:${RUN_DIR}" \
$mount_clause \
--cap-add SYS_ADMIN \
--device /dev/fuse \
$RUFUS_DOCKER_IMAGE \
tail -f /dev/null) # TODO: do I need this line

EOF
    fi
}
export -f write_container_start_piece

write_bwa_index_pieces() {
    local out_script="$1"

    exist=$(check_for_bwa_indexes)
    if [ "${exist}" == "false" ]; then

        cat <<EOF >> "$out_script"
# Generate BWA indexes for reference fasta
# Note: this step can be skipped next time you run by copying indexes into same location of $REFERENCE_FASTA
EOF

        if [ "$container_type" == "$SINGULARITY" ]; then

        # TODO: what is in-dir and out here? should already have dirs bound from start?
            cat <<EOF >> "$out_script" 
srun --nodes=1 --ntasks=1 --nodelist=\${FIRST_NODE}" \
     singularity exec instance://rufus_instance_\${SLURM_NODEID} bash -lc "echo 'Building reference indexes on ' \$(hostname); ${INDEX_BUILD_SCRIPT} --in-dir /shared/path --out /shared/path/final_output"
EOF
        else
            echo -e "docker exec \${CONTAINER_ID} bash ${INDEX_BUILD_SCRIPT} ${REFERENCE_FASTA}" >> "$out_script"
        fi
    fi
}
export -f write_bwa_indexes

write_rufus_execution_piece() {
    local out_script="$1"

    # Some args will be region agnostic
    rufus_args=$(get_rufus_args)

    if [ "$WINDOWED_MODE" ]; then

        cat <<EOF >> "$out_script"
# Get regional specific run-time args
region=$(head -n "\$(\$SLURM_NODEID + 1)" "$REGION_FILE" | tail -n 1)
hash_arg=\$(singularity exec instance://rufus_instance_\${SLURM_NODEID} "$HASH_SCRIPT" "\$region" "kg1")

EOF
        # Get spec args
        spec_args=$(get_slurm_specs "$WINDOWED_MODE" "$container_type")

        if [ "$container_type" == "$SINGULARITY" ]; then
            # Regional, using singularity
            ntasks="${spec_args[0]}"
            cpus_per_task="${spec_args[1]}"
            mem_per_task="${spec_args[2]}"

            cat <<EOF >> "$out_script"
srun --ntasks=${ntasks} \
     --cpus-per-task=${cpus_per_task} \
     --mem=${mem_per_task} \
     --kill-on-bad-exit=1 \
     --output=task_logs/task_%t_%N_%j.out \
     bash -lc "singularity exec instance://rufus_instance_\${SLURM_NODEID} $ENTRY_SCRIPT $rufus_args \$hash_arg -r \$region"
EOF

        else
            job_phrase=""
            if [ -n "$PARALLEL_JOBS" ]; then
                job_phrase="-j $PARALLEL_JOBS"
            fi

            # Regional, using docker
            # TODO: how should I access PR_WRAPPER? Pull it on to host machine?
            cat <<EOF >> "$out_script"
parallel "$job_phrase" bash "${PR_WRAPPER}" "\$CONTAINER_ID" {} :::: "$REGION_FILE"
EOF
        fi
    else 
        cat <<EOF >> "$out_script"
RUFUS_CMD="$ENTRY_SCRIPT \
$rufus_args"

EOF
        if [ "$container_type" == "$SINGULARITY" ]; then
            # Full genome, no parallelism, using singularity
            cat <<EOF >> "$out_script"
srun singularity exec "instance://rufus_instance_\${SLURM_NODEID}" "bash -lc \
    \$RUFUS_CMD"
EOF
        else
            # Full genome, no parallelism, using docker
            cat <<EOF >> "$out_script"
docker exec "\$CONTAINER_ID" "bash -lc \
    \$RUFUS_CMD \
    > ${LOGS_DIR}wg.out \
    2> ${LOGS_DIR}wg.err"
EOF
        fi
    fi
}
export -f write_rufus_execution_piece

write_post_process_piece() {
    local out_script="$1"
    local post_args=""
    post_args=$(get_post_string)

    if [ "$container_type" == "$SINGULARITY" ]; then
        # TODO: verify that I don't need to add a --depend arg here on the main parallel execution

        cat <<EOF >> "$out_script"
FIRST_NODE=\$(scontrol show hostnames "\$SLURM_JOB_NODELIST" | head -n1)
echo "Running final aggregation on first node: \${FIRST_NODE}"

srun --nodes=1 --ntasks=1 --nodelist="\${FIRST_NODE}" \
     singularity exec instance://rufus_instance_\${SLURM_NODEID} bash -lc "echo 'Final aggregation running on ' \$(hostname); ${POST_SCRIPT} ${post_args} --in-dir /shared/path --out /shared/path/final_output"

EOF
    else
        echo -e "docker exec \${CONTAINER_ID} bash ${POST_SCRIPT} ${post_args}" >> "$out_script"
    fi
}
export -f write_post_process_piece

write_container_stop_piece() {
    local out_script="$1"

    if [ "$container_type" == "$SINGULARITY" ]; then

        cat <<EOF >> "$out_script"
echo "Stopping instances on all nodes..."
srun --ntasks-per-node=1 --exclusive \
     bash -lc "singularity instance stop "rufus_instance_\${SLURM_NODEID}" || echo 'stop failed or already stopped on' \$(hostname)"

EOF
    else

        cat <<EOF >> "$out_script"
echo "Stopping docker container..."
echo docker stop "${CONTAINER_NAME}" >/dev/null

EOF
    fi
}
export -f write_container_stop_piece

# Writes entire slurm or bash script that user will execute to launch RUFUS
write_out_script() {
    if [ "$container_type" == "$SINGULARITY" ]; then
        out_script="${WORKING_DIR}/run_rufus_singularity.slurm"
    else
        out_script="${WORKING_DIR}/run_rufus_docker.sh"
    fi
        write_header "$out_script"
        write_container_start_piece "$out_script"

        cat <<EOF >> "$out_script"
start_time=\$(date +%s)" 
"absolute_start=\$start_time"

EOF
        write_bwa_index_piece "$out_script"

        cat <<EOF >> "$out_script"
end_time=\$(date +%s)
elapsed=\$((end_time - start_time))
printf -v human "%02d:%02d:%02d" \$((elapsed/3600)) \$(((elapsed%3600)/60)) \$((elapsed%60))
echo "BWA index creation complete. Step run time: \$human"
start_time=\$(date +%s)

EOF

        write_rufus_execution_piece "$out_script"

        cat <<EOF >> "$out_script"
end_time=\$(date +%s)
elapsed=\$((end_time - start_time))
printf -v human "%02d:%02d:%02d" \$((elapsed/3600)) \$(((elapsed%3600)/60)) \$((elapsed%60))
echo "RUFUS calling stage complete. Step run time: \$human"
start_time=\$(date +%s)

EOF
        write_post_process_piece "$out_script"

        cat <<EOF >> "$out_script"
end_time=\$(date +%s)
elapsed=\$((end_time - start_time))
printf -v human "%02d:%02d:%02d" \$((elapsed/3600)) \$(((elapsed%3600)/60)) \$((elapsed%60))
echo "RUFUS post-processing stage complete. Step run time: \$human"

EOF

        write_container_stop_piece "$out_script"

        cat <<EOF >> "$out_script"
end_time=\$(date +%s)
elapsed=\$((end_time - absolute_start))
printf -v human "%02d:%02d:%02d" \$((elapsed/3600)) \$(((elapsed%3600)/60)) \$((elapsed%60))
Entire RUFUS run time: \$human
EOF
}
export -f write_out_script


#------------------- DO WORK -------------------#

# Parse config file and write to temp.env file
parse_config_file

# Write out launch script
write_out_script