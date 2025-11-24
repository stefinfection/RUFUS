# Generates launch script that using either GNU Parallel or Singularity

container_type="$1"

# Constants
config_file=/temp/rufus_config.yaml
SINGULARITY="singularity"
DOCKER="docker"

# Parse config file and write to temp.env file
parse_config_file

exist=$(check_for_bwa_indexes)
if [ ${exist} == "false" ]; then
    echo "# Generate BWA indexes for reference fasta" >> "$launch_out"
    echo "# Note: this step can be skipped next time you run by copying indexes into same location of $REFERENCE_FASTA" >> "$launch_out"
    echo "docker exec ${CONTAINER_ID} bash /opt/RUFUS/resource_helpers/build_bwa_indexes.sh ${REFERENCE_FASTA}" >> "$launch_out"
fi

# Returns mount clause for controls, 1000G, subject, and references/indexes
get_container_mount_clause() {
    # Get option flag based on container type
    mount_opt="-v"
    if [ "$container_type" == "$SINGULARITY" ]; then
        mount_opt="--bind"
    fi

    mount_clause=""

    # Add subject file + index
    subject_mount_array=$(set_up_subject) # TODO: Return both cram/bam AND crai/bai here
    for sub in "${subject_mount_array[@]}"; do
        sub_basename=$(basename $sub)
        mount_clause+="$mount_opt $sub:/mnt/rufus_temp/$sub_basename:ro"
    done

    control_mount_array=$(set_up_controls) # TODO: Return both cram/bam AND crai/bai here
    for ctrl in "${control_mount_array[@]}"; do
        ctrl_basename=$(basename $ctrl)
        mount_clause+="$mount_opt $ctrl:/mnt/rufus_temp/$ctrl_basename:ro "
    done
    
    kg1_mount=$(set_up_kg1)
    mount_clause+="$mount_opt $kg1_mount"

    ref_mount_array=$(set_up_ref) # TODO: return all indexes too if they exist
    for ref in "${ref_mount_array[@]}"; do
        mount_clause+="$mount_opt $ref"
    done

    echo "$mount_clause"
}

# TODO: should this be write out main portion and then have a nother function to write out post-process piece?
write_out_script() {
    if [ "$container_type" == "$SINGULARITY" ]; then
        out_script="run_rufus_singularity.slurm"

        # get #SBATCH header file for job run script (TODO: how to code this human readable that will echo correctly)

        # write out regional specific arg fxn calls
        # It gets the regional specific arguments for control and kg1 hashes (these functions will be internal to rufus container now)
        # It sruns the individual region jobs
        
        # write post-process script
        # slurm_id=sbatch slurm script
        # sbatch post-process script --depend-on:slurm_id



    else
        out_script="run_rufus_docker.sh"

        echo "$CONTAINER_ID=$(docker run -d --rm --name rufus-worker \
            -u "$(id -u):$(id -g)" \
            -v "${WORKING_DIR}:/mnt" \
            --cap-add SYS_ADMIN \
            --device /dev/fuse \" >> $output_script
        
        local mount_clause=$(get_container_mount_clause)
        echo "$mount_clause"
        
        echo "$RUFUS_DOCKER_IMAGE" \
            tail -f /dev/null)" >> $output_script
    fi
        echo "start_time=$(date +%s)" >> $out_script

}

# Write out bash script
    # get_container_mount_clause
    # get_input_files
    # start docker container
    # pull out region worker/internal launch helpers script
    # write out arguments to rufus.cmd
    # parallel run - might still need a wrapper here to coordinate proper error/out files, but should be able to use same get/fetch hash fxns
    # docker exec post process
    # stop docker container
    # Output total time taken


# Returns the clause to add after `singularity instance start`
get_container_start_args(cont_type) {

}

get_slurm_script () {
       

}
# This function outputs the script that will be launched once on each node
# Divide region input file by number of nodes to get processing regions per node - just get range relative to entire file so don't have to reproduce
# Make const sbatch header + fill in num nodes + slurm info
# Make sure to use --nodes=NUM_NODES according to env file instead of arrays
# get start command for singularity get_container_start_args(singularity)
# start container with srun
# write out arguments to rufus.cmd
# use srun to start array jobs, instead of writing another script to sbatch (just fill in cpus-per-task and mem to have slurm max parallelism)
# try 8G mem and 10 cpus per task (so 20 threads to rufus) if region - otherwise do whole genome
# If first script in all nodes queued (FIRST_NODE=$(scontrol show hostnames "$SLURM_JOB_NODELIST" | head -n1)), then run post-process

# If attempt to use more nodes than region (in case doing whole genome mode), warn only using one

# get control argument formatted for post process (will use for both so do before if statement)