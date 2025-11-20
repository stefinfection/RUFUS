# if building refs, start slurm job for this on a single node
# set flag that start scripts all need to wait on this job

# if singularity:
    # get_slurm_script
    # write post-process script
    # slurm_id=sbatch slurm script
    # sbatch post-process script --depend-on:slurm_id