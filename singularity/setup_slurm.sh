#!/bin/bash
# This script creates three files in the directory it's run in (pwd):
# 1. A rufus call slurm script
# 2. A rufus post-process slurm script
# 3. A bash script to batch submit the two above slurm scripts

#LOCAL_TESTING_UTIL_PATH=/home/ubuntu/RUFUS/singularity/launch_utilities/
#UTIL_PATH=$LOCAL_TESTING_UTIL_PATH
UTIL_PATH=/opt/RUFUS/singularity/launch_utilities/

PARSER=${UTIL_PATH}arg_parser.sh
. $PARSER "$@"

GENOME_HELPERS=${UTIL_PATH}genome_helpers.sh
. $GENOME_HELPERS

CHUNK_UTILITIES=${UTIL_PATH}chunk_utilities.sh
. $CHUNK_UTILITIES

NUM_CHUNKS=$(get_num_chunks $WINDOW_SIZE_RUFUS_ARG $GENOME_BUILD_RUFUS_ARG)

WORKING_DIR=$(pwd)

echo -en "##RUFUS_callCommand=" > rufus.cmd

# TODO: if ARRAY_SIZE_LIMIT < NUM_CHUNKS need to split these into multiple batch scripts
# Shift to 0-based slurm array index
# ADJ_SLURM_ARRAY_JOB_LIMIT_RUFUS_ARG=$(($SLURM_ARRAY_JOB_LIMIT_RUFUS_ARG - 1))
#if [ "$ADJ_SLURM_ARRAY_JOB_LIMIT_RUFUS_ARG" -lt "$NUM_CHUNKS" ]; then

	# TODO: left off here
	# TODO: get number of iterations by ceil(NUM_CHUNKS/ADJ_LIMIG)		
	# TODO: need to wrap below script composing in for loop and label each script by index number
	# TODO: also shift NUM_CHUNKS by one if not (to include 0)
#fi

# Compose run script(s)
RUFUS_SLURM_SCRIPT="rufus_call.slurm"
HEADER_LINES=("#!/bin/bash"
"#SBATCH --job-name=rufus_call"
"#SBATCH --time=${SLURM_TIME_LIMIT_RUFUS_ARG}" 
"#SBATCH --account=${SLURM_ACCOUNT_RUFUS_ARG}" 
"#SBATCH --partition=${SLURM_PARTITION_RUFUS_ARG}"
"#SBATCH --cpus-per-task=${THREAD_LIMIT_RUFUS_ARG}"
"#SBATCH -o ${WORKING_DIR}/slurm_out/%A_%a.out"
"#SBATCH -e ${WORKING_DIR}/slurm_err/%A_%a.err"
)

# Don't overwrite a run if already exists
if [ -f "$RUFUS_SLURM_SCRIPT" ]; then
	echo "ERROR: $RUFUS_SLURM_SCRIPT already exists - are you overwriting an existing output? Please delete run_rufus.slurm and retry"
	exit 1
fi

for line in "${HEADER_LINES[@]}"
do
    echo -e "$line" >> $RUFUS_SLURM_SCRIPT
done

if [ ! -z $EMAIL_RUFUS_ARG ]; then
	echo -e "#SBATCH --mail-type=ALL" >> $RUFUS_SLURM_SCRIPT
	echo -e "#SBATCH --mail-user=${EMAIL_RUFUS_ARG}" >> $RUFUS_SLURM_SCRIPT
fi

if [ "$WINDOW_SIZE_RUFUS_ARG" = "0" ]; then
	echo "" >> $RUFUS_SLURM_SCRIPT
	echo -e "REGION_ARG=\"\"" >> $RUFUS_SLURM_SCRIPT
else
	echo -e "#SBATCH -a 0-${NUM_CHUNKS}%${SLURM_JOB_LIMIT_RUFUS_ARG}" >> $RUFUS_SLURM_SCRIPT
	echo "" >> $RUFUS_SLURM_SCRIPT
	echo -e "region_arg=\$(singularity exec ${CONTAINER_PATH_RUFUS_ARG}/rufus.sif bash /opt/RUFUS/singularity/launch_utilities/get_region.sh \"\$SLURM_ARRAY_TASK_ID\" \"$WINDOW_SIZE_RUFUS_ARG\" \"$GENOME_BUILD_RUFUS_ARG\")" >> $RUFUS_SLURM_SCRIPT
	echo -e "REGION_ARG=\"-R \$region_arg\"" >> $RUFUS_SLURM_SCRIPT
fi

echo -en "srun --mem=0 singularity exec --bind ${HOST_DATA_DIR_RUFUS_ARG}:/mnt ${CONTAINER_PATH_RUFUS_ARG}/rufus.sif bash /opt/RUFUS/runRufus.sh -s /mnt/$SUBJECT_RUFUS_ARG " >> $RUFUS_SLURM_SCRIPT
echo -en "srun --mem=0 singularity exec --bind ${HOST_DATA_DIR_RUFUS_ARG}:/mnt ${CONTAINER_PATH_RUFUS_ARG}/rufus.sif bash /opt/RUFUS/runRufus.sh -s /mnt/$SUBJECT_RUFUS_ARG " >> rufus.cmd
for control in "${CONTROLS_RUFUS_ARG[@]}"; do
	echo -en "-c $control " >> $RUFUS_SLURM_SCRIPT
	echo -en "-c $control " >> rufus.cmd
done

if [ ! -z "$REF_HASH_RUFUS_ARG" ]; then
	echo -en "-f $REF_HASH_RUFUS_ARG " >> $RUFUS_SLURM_SCRIPT
	echo -en "-f $REF_HASH_RUFUS_ARG " >> rufus.cmd
fi
echo -e "-r $REFERENCE_RUFUS_ARG -m $KMER_DEPTH_CUTOFF_RUFUS_ARG -k 25 -t $THREAD_LIMIT_RUFUS_ARG -L -vs \$REGION_ARG" >> $RUFUS_SLURM_SCRIPT
echo -e "-r $REFERENCE_RUFUS_ARG -m $KMER_DEPTH_CUTOFF_RUFUS_ARG -k 25 -t $THREAD_LIMIT_RUFUS_ARG -L -vs \$REGION_ARG" >> rufus.cmd

# Compose post-process slurm wrapper
PP_SLURM_SCRIPT="rufus_post_process.slurm" # Slurm wrapper for post process script
PP_HEADER_LINES=("#!/bin/bash"
"#SBATCH --job-name=rufus_post_process"   
"#SBATCH --account=${SLURM_ACCOUNT_RUFUS_ARG}" 
"#SBATCH --partition=${SLURM_PARTITION_RUFUS_ARG}"
"#SBATCH --output=${WORKING_DIR}/slurm_out/rufus_post_process_%j.out"   
"#SBATCH --error=${WORKING_DIR}/slurm_err/rufus_post_process_%j.err" 
"#SBATCH --nodes=1"
)

for line in "${PP_HEADER_LINES[@]}"
do
    echo -e "$line" >> $PP_SLURM_SCRIPT
done

if [ ! -z $EMAIL_RUFUS_ARG ]; then
    echo -e "#SBATCH --mail-type=ALL" >> $PP_SLURM_SCRIPT
    echo -e "#SBATCH --mail-user=${EMAIL_RUFUS_ARG}" >> $PP_SLURM_SCRIPT
fi
echo "" >> $PP_SLURM_SCRIPT

IFS=$','
CONTROL_STRING="${CONTROLS_RUFUS_ARG[*]}"

BOUND_DATA_DIR="/mnt"

echo -e "srun singularity exec --bind ${HOST_DATA_DIR_RUFUS_ARG}:/mnt ${CONTAINER_PATH_RUFUS_ARG}/rufus.sif bash /opt/RUFUS/post_process/post_process.sh -w $WINDOW_SIZE_RUFUS_ARG -r $REFERENCE_RUFUS_ARG -c $CONTROL_STRING -s $SUBJECT_RUFUS_ARG -d ${BOUND_DATA_DIR}" >> $PP_SLURM_SCRIPT

echo -en "##RUFUS_postProcessCommand=" >> rufus.cmd
echo -e "srun singularity exec --bind ${HOST_DATA_DIR_RUFUS_ARG}:/mnt ${CONTAINER_PATH_RUFUS_ARG}/rufus.sif bash /opt/RUFUS/post_process/post_process.sh -w $WINDOW_SIZE_RUFUS_ARG -r $REFERENCE_RUFUS_ARG -c $CONTROL_STRING -s $SUBJECT_RUFUS_ARG -d ${BOUND_DATA_DIR}" >> rufus.cmd
mv rufus.cmd ${HOST_DATA_DIR_RUFUS_ARG}

# Compose invocation script to be executed outside of container
EXE_SCRIPT=launch_rufus.sh
echo -e "#!/bin/bash" > $EXE_SCRIPT
echo -e "" >> $EXE_SCRIPT
echo -e "# This script should be executed after calling the container setup_slurm.sh helper. It requires $PP_SLURM_SCRIPT and $RUFUS_SLURM_SCRIPT to be present in the same directory." >> $EXE_SCRIPT 
echo -e "# Insert command for your system to load singularity here if needed (e.g. module load singularity)"
echo "" >> $EXE_SCRIPT
echo -e "# Launch calling job" >> $EXE_SCRIPT
echo -e "ARRAY_JOB_ID=\$(sbatch --parsable $RUFUS_SLURM_SCRIPT)" >> $EXE_SCRIPT
echo -e "" >> $EXE_SCRIPT
echo -e "# Launch post-process job - will wait on calling phase to complete" >> $EXE_SCRIPT
echo -e "sbatch --depend=afterany:\$ARRAY_JOB_ID $PP_SLURM_SCRIPT" >> $EXE_SCRIPT

echo -e "Slurm scripts ready to execute with $EXE_SCRIPT. Please make sure singularity is available in your environment, and then run... "
echo -e "bash $EXE_SCRIPT"
