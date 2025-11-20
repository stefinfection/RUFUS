# Generates launch script that uses GNU parallel; runs inside container

config_file=/temp/rufus_config.yaml

# Parse config file
# Import bash helper module

# check_for_bwa_indexes()
    # TODO: when do I call this - do I just add it to the bash script at the top?
    # notify user that index refs are being built and that this step can be skipped next time to save time

# get_container_mount_clause()
    # IFS='|' read -r ref_mount build_refs < <(set_up_ref)
    # control_mount=$(set_up_controls) || exit 1
    # kg1_mount=$(set_up_kg1) || exit 1 
    # input_mount_clause="$ref_mount $control_mount $kg1_mount -v ${subject_path}:/mnt/rufus_temp/${subject_base}:ro "

# get_input_files (subject + check for index)
    # subject_path=$(realpath "${SUBJECT_FILE}")
    # subject_base=$(basename "${subject_path}")

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
get_container_start_args(cont_type) {}


# It starts the singularity container
# It gets the regional specific arguments for control and kg1 hashes (these functions will be internal to rufus container now)
# It sruns the individual region jobs
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