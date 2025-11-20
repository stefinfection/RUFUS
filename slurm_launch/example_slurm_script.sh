#!/bin/bash
#SBATCH --job-name=sing_inst_4nodes
#SBATCH --nodes=4                    # change to desired node count
#SBATCH --time=1-00:00:00
#SBATCH --partition=normal
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err

### -------- USER VARIABLES: edit these --------
IMAGE="/path/to/my_image.sif"               # shared container image path
INSTANCE_NAME="my_instance"                 # name used for the instance on each node
COMMANDS_FILE="/shared/path/commands.txt"   # shared file: one command per line
FINAL_SCRIPT="/shared/path/final_aggregate.sh" # final aggregation script (must be on shared FS)
CPUS_PER_TASK=1                              # CPU cores per small task
MEM_PER_TASK=1G                              # memory per small task
### --------------------------------------------

# todo: adapt to regions file instead of commands file

# derived
TOTAL_TASKS=$(wc -l < "${COMMANDS_FILE}")
echo "Job started on $SLURM_JOB_NUM_NODES nodes (requested). Total tasks: ${TOTAL_TASKS}"

# Prepare directories
mkdir -p logs task_logs
echo "Logs => $(realpath logs) ; task logs => $(realpath task_logs)"

# Safety: trap to stop instances on exit (best-effort)
cleanup() {
  echo "Cleanup: stopping instances on all nodes..."
  # attempt to stop on every allocated node
  srun --ntasks-per-node=1 --exclusive \
       bash -lc "singularity instance stop ${INSTANCE_NAME} || true"
}
trap cleanup EXIT

# 1) Start one instance per node (blocking: returns after instance started on each node)
echo "Starting one instance per node..."
srun --ntasks-per-node=1 --exclusive \
     bash -lc "hostname; singularity instance start --writable-tmpfs ${IMAGE} ${INSTANCE_NAME} && echo 'instance started' || (echo 'instance START FAILED on' \$(hostname); exit 1)"

# quick verification (optional)
echo "Verifying instance on each node..."
srun --ntasks-per-node=1 --exclusive \
     bash -lc "echo 'on' \$(hostname):; singularity instance list | sed -n '1,100p' ; singularity exec instance://${INSTANCE_NAME} hostname || echo 'exec failed on' \$(hostname)"

# 2) Create a per-task worker script in the job's temp dir to avoid nested-quoting issues
WORKER_SH="${SLURM_TMPDIR:-/tmp}/worker_${SLURM_JOB_ID}.sh"
cat > "${WORKER_SH}" <<'__WORKER_SH__'
#!/bin/bash
set -euo pipefail

# this worker script is executed under srun for each task
# SLURM provides SLURM_PROCID (0..N-1) for the task id in the job-step
if [ -z "${SLURM_PROCID+x}" ]; then
  echo "SLURM_PROCID not set. Are you running under srun?" >&2
  exit 2
fi

IDX=$(( SLURM_PROCID + 1 ))   # line number (1-based)
CMD="$(sed -n "${IDX}p" "${COMMANDS_FILE}")"

echo "[$(date)] [node=$(hostname)] procid=${SLURM_PROCID} idx=${IDX} running: ${CMD}"
# run inside the local instance on this node
singularity exec instance://${INSTANCE_NAME} bash -lc "${CMD}"
RET=$?
echo "[$(date)] [node=$(hostname)] procid=${SLURM_PROCID} idx=${IDX} exit=${RET}"
exit ${RET}
__WORKER_SH__

chmod +x "${WORKER_SH}"
echo "Worker script written to ${WORKER_SH}"

# export vars so the worker script can read them (srun will propagate environment)
export COMMANDS_FILE INSTANCE_NAME

# 3) Launch all small tasks with a single blocking srun that spans the allocation
#    This srun will not return until every task on every node completes.
echo "Launching ${TOTAL_TASKS} small tasks across the allocation (this will block until all finish)..."
srun --ntasks=${TOTAL_TASKS} \
     --cpus-per-task=${CPUS_PER_TASK} \
     --mem=${MEM_PER_TASK} \
     --kill-on-bad-exit=1 \
     --output=task_logs/task_%t_%N_%j.out \
     "${WORKER_SH}"

echo "All small tasks finished across all nodes."

# 4) Run final aggregation ON THE FIRST ALLOCATED NODE
FIRST_NODE=$(scontrol show hostnames "$SLURM_JOB_NODELIST" | head -n1)
echo "Running final aggregation on first node: ${FIRST_NODE}"

srun --nodes=1 --ntasks=1 --nodelist="${FIRST_NODE}" \
     singularity exec instance://${INSTANCE_NAME} bash -lc "echo 'Final aggregation running on ' \$(hostname); ${FINAL_SCRIPT} --in-dir /shared/path --out /shared/path/final_output"

echo "Final aggregation finished."

# 5) Stop instances on all nodes (cleanup will also run on EXIT)
echo "Stopping instances on all nodes..."
srun --ntasks-per-node=1 --exclusive \
     bash -lc "singularity instance stop ${INSTANCE_NAME} || echo 'stop failed or already stopped on' \$(hostname)"

# disable trap cleanup (we already cleaned up)
trap - EXIT
echo "Done. Instances stopped. Job complete."
