#!/bin/bash
CONTAINER_PATH=$1

PARSER=${CONTAINER_PATH}/opt/RUFUS/singularity/launch_utilities/arg_parser.sh
. $PARSER "$@"

GENOME_HELPERS=/opt/RUFUS/singularity/launch_utilities/genome_helpers.sh
. $GENOME_HELPERS

CHUNK_UTILITIES=/opt/RUFUS/singularity/launch_utilites/chunk_utilites.sh
. $CHUNK_UTILITIES
NUM_CHUNKS=$(get_num_chunks.sh $WINDOW_SIZE_RUFUS_ARG $GENOME_BUILD_RUFUS_ARG)

WORKING_DIR=$(pwd)
mkdir slurm_out

if [ "$WINDOW_SIZE_RUFUS_ARG" = "0" ]; then
	mkdir -p slurm_out/individual_jobs
fi

# Compose run script
RUFUS_SLURM_SCRIPT="run_rufus.slurm"
HEADER_LINES=("#!/bin/bash"
"#SBATCH --job-name=RUFUS"
"#SBATCH --time=${SLURM_TIME_LIMIT_RUFUS_ARG}" 
"#SBATCH --account=${SLURM_ACCOUNT_RUFUS_ARG}" 
"#SBATCH --partition=${SLURM_PARTITION_RUFUS_ARG}"
"#SBATCH --nodes=1"
"#SBATCH -o ${WORKING_DIR}/slurm_out/%j.out"
"#SBATCH -e ${WORKING_DIR}/slurm_err/%j.err"
"#SBATCH --mail-type=ALL" 
"#SBATCH --mail-user=${EMAIL_RUFUS_ARG}"
"#SBATCH -a 0-${NUM_CHUNKS}%${SLURM_JOB_LIMIT_RUFUS_ARG}"
)

# Don't overwrite a run if already exists
if [ -f "$RUFUS_SLURM_SCRIPT" ]; then
	echo "ERROR: run_rufus.slurm already exists - are you overwriting an existing output? Please delete run_rufus.slurm and retry"
	exit 1
fi

for line in "${HEADER_LINES[@]}"
do
    echo -e "$line" >> $RUFUS_SLURM_SCRIPT
done
echo "\n" >> $RUFUS_SLURM_SCRIPT

# todo: need to test getChunk.sh
echo "REGION_ARG=$(bash getChunk.sh "${SLURM_ARRAY_TASK_ID}")" >> $RUFUS_SLURM_SCRIPT

echo -en "singularity exec --bind ${HOST_DATA_DIR_RUFUS_ARG}:/mnt ${CONTAINER_PATH_RUFUS_ARG}/rufus.sif bash /opt/RUFUS/runRufus.sh -s /mnt/$SUBJECT_RUFUS_ARG " >> $RUFUS_SLURM_SCRIPT
for control in "${CONTROLS_RUFUS_ARG[@]}"; do
	echo -en "-c $control "
done

if [ -z "$REF_HASH_RUFUS_ARG" ]; then
	echo -en "-f $REF_HASH_RUFUS_ARG " >> $RUFUS_SLURM_SCRIPT
fi
# todo: I won't have slurm_array_task_id here yet...
echo -e "-r $REFERENCE_RUFUS_ARG -m $KMER_DEPTH_CUTOFF_RUFUS_ARG -k 25 -t $THREAD_LIMIT_RUFUS_ARG -L -vs $REGION_ARG -v $SLURM_ARRAY_TASK_ID" >> $RUFUS_SLURM_SCRIPT

echo "testing outs:"
echo "$CONTROL_STRING_RUFUS_ARG"
SUBJECT_STRING=$(basename $SUBJECT_RUFUS_ARG)
echo "$SUBJECT_STRING"

# EXECUTE BASH SCRIPT
#sbatch $RUFUS_SLURM_SCRIPT

#todo: what args do I need to pass here?
#SUBJECT_STRING=$(basename $SUBJECT_RUFUS_ARG)
#sbatch --depend=afterany:$SLURM_ARRAY_JOB_ID singularity exec bash postProcessRufus.sh -w $WINDOW_SIZE_RUFUS_ARG -r $REFERENCE_RUFUS_ARG -c $CONTROL_STRING_RUFUS_ARG -s $SUBJECT_STRING -d $WORKING_DIR 
