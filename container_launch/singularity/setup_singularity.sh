# if building refs, start slurm job for this on a single node
# set flag that start scripts all need to wait on this job

# if singularity:
    # get_slurm_script
    # write post-process script
    # slurm_id=sbatch slurm script
    # sbatch post-process script --depend-on:slurm_id

    # It starts the singularity container
# It gets the regional specific arguments for control and kg1 hashes (these functions will be internal to rufus container now)
# It sruns the individual region jobs