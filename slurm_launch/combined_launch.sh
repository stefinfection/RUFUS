#!/bin/bash

# import launch helpers OR put all of them in here so we don't have to have user download another file

# Check for env file

# Set working directory

# Only helper functions outside of container are those necessary to START the container, not exec individual jobs
# container_type=check_req_variables()
# IFS='|' read -r ref_mount build_refs < <(set_up_ref)
# notify user that index refs are being built and that this step can be skipped next time to save time
# if building refs, start slurm job for this on a single node
# set flag that start scripts all need to wait on this job
# TODO: do I start container and do this? have to if we want to send it to slurm

# control_mount=$(set_up_controls) || exit 1
# kg1_mount=$(set_up_kg1) || exit 1
# subject_path=$(realpath "${SUBJECT_FILE}")
# subject_base=$(basename "${subject_path}")
# input_mount_clause="$ref_mount $control_mount $kg1_mount -v ${subject_path}:/mnt/rufus_temp/${subject_base}:ro "

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

# if singularity:
    # get_slurm_script
    # write post-process script
    # slurm_id=sbatch slurm script
    # sbatch post-process script --depend-on:slurm_id

#if docker:
    # get_docker_start_args
    # start docker container
    # pull out region worker/internal launch helpers script
    # write out arguments to rufus.cmd
    # parallel run - might still need a wrapper here to coordinate proper error/out files, but should be able to use same get/fetch hash fxns
    # docker exec post process
    # stop docker container

# Output total time taken

