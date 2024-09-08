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
if [ -f "rufus_call*.slurm" ]; then
    echo "ERROR: a rufus_call.slurm script already exists - are you overwriting an existing output? Please delete run_rufus*.slurm and retry"
    exit 1
fi

# Compose run script(s)
echo -en "##RUFUS_callCommand=" > rufus.cmd

# Figure out how many call scripts we need
NUM_CALL_SLURMS=$(( NUM_CHUNKS / SLURM_ARRAY_JOB_LIMIT_RUFUS_ARG ))
REM=$(( NUM_CHUNKS % SLURM_ARRAY_JOB_LIMIT_RUFUS_ARG ))
if [ "$REM" -gt 0 ]; then
    NUM_CALL_SLURMS=$(( NUM_CALL_SLURMS + 1 ))
fi

# If doing non-windowed run, or fits into array, similar composure
SINGLE_MODE="FALSE"
CALL_SLURM_SCRIPT="rufus_scaffold.slurm"
if [ "$WINDOW_SIZE_RUFUS_ARG" = "0" ] || [ "$NUM_CALL_SLURMS" = "1" ]; then
	SINGLE_MODE="TRUE"
	CALL_SLURM_SCRIPT="rufus_call.slurm"
fi	

# Compose the call script
for line in "${HEADER_LINES[@]}"
do
   	echo -e "$line" >> $CALL_SLURM_SCRIPT
done

if ["$SINGLE_MODE" = "FALSE" ];
	echo -e "REGION_INDEX=\$1" >> $CALL_SLURM_SCRIPT
fi

if [ ! -z $EMAIL_RUFUS_ARG ]; then
   	echo -e "#SBATCH --mail-type=ALL" >> $CALL_SLURM_SCRIPT
   	echo -e "#SBATCH --mail-user=${EMAIL_RUFUS_ARG}" >> $CALL_SLURM_SCRIPT
fi		

if [ "$WINDOW_SIZE_RUFUS_ARG" = "0" ]; then
   	echo "" >> $CALL_SLURM_SCRIPT
   	echo -e "REGION_ARG=\"\"" >> $CALL_SLURM_SCRIPT
else
	ADJUSTED_END=$(( SLURM_ARRAY_JOB_LIMIT_RUFUS_ARG - 1 ))
	
	# Only use arrays if we fit into a single script
	if [ "$SINGLE_MODE" = "TRUE" ]; then
   		echo -e "#SBATCH -a 0-${ADJUSTED_END}%${SLURM_JOB_LIMIT_RUFUS_ARG}" >> $CALL_SLURM_SCRIPT
   		echo "" >> $CURR_CALL_SCRIPT
   		echo -e "region_arg=\$(singularity exec ${CONTAINER_PATH_RUFUS_ARG} bash /opt/RUFUS/singularity/launch_utilities/get_region.sh \"\$SLURM_ARRAY_TASK_ID\" \"$WINDOW_SIZE_RUFUS_ARG\" \"$GENOME_BUILD_RUFUS_ARG\")" >> $CURR_CALL_SCRIPT
	else
		echo "" >> $CURR_CALL_SCRIPT
        echo -e "region_arg=\$(singularity exec ${CONTAINER_PATH_RUFUS_ARG} bash /opt/RUFUS/singularity/launch_utilities/get_region.sh \"\$REGION_INDEX\" \"$WINDOW_SIZE_RUFUS_ARG\" \"$GENOME_BUILD_RUFUS_ARG\")" >> $CURR_CALL_SCRIPT
	fi

   	echo -e "REGION_ARG=\"-R \$region_arg\"" >> $CALL_SLURM_SCRIPT
fi

echo -en "srun --mem=0 singularity exec --bind ${HOST_DATA_DIR_RUFUS_ARG}:/mnt ${CONTAINER_PATH_RUFUS_ARG} bash /opt/RUFUS/runRufus.sh -s /mnt/$SUBJECT_RUFUS_ARG " >> $CALL_SLURM_SCRIPT
echo -en "srun --mem=0 singularity exec --bind ${HOST_DATA_DIR_RUFUS_ARG}:/mnt ${CONTAINER_PATH_RUFUS_ARG} bash /opt/RUFUS/runRufus.sh -s /mnt/$SUBJECT_RUFUS_ARG " >> rufus.cmd

# Add control args
for control in "${CONTROLS_RUFUS_ARG[@]}"; do
	echo -en "-c $control " >> $CALL_SLURM_SCRIPT
   	echo -en "-c $control " >> rufus.cmd
done	

# Add in optional hashes if provided
if [ ! -z "$REFERENCE_HASH_RUFUS_ARG" ]; then
  	echo -en "-f $REFERENCE_HASH_RUFUS_ARG " >> $CALL_SLURM_SCRIPT
   	echo -en "-f $REFERENCE_HASH_RUFUS_ARG " >> rufus.cmd
fi
if [ ! -z "$EXCLUDE_HASH_LIST_RUFUS_ARG" ]; then
   	for exclude in "${EXCLUDE_HASH_LIST_RUFUS_ARG[@]}"; do
    	echo -en "-e $exclude " >> $CALL_SLURM_SCRIPT
       	echo -en "-e $exclude " >> rufus.cmd
   	done
fi

echo -e "-r $REFERENCE_RUFUS_ARG -m $KMER_DEPTH_CUTOFF_RUFUS_ARG -k 25 -t $THREAD_LIMIT_RUFUS_ARG -L -vs \$REGION_ARG" >> $CALL_SLURM_SCRIPT
echo -e "-r $REFERENCE_RUFUS_ARG -m $KMER_DEPTH_CUTOFF_RUFUS_ARG -k 25 -t $THREAD_LIMIT_RUFUS_ARG -L -vs \$REGION_ARG" >> rufus.cmd

# Make throttled slurm script to run scaffold	
if [ "$SINGLE_MODE" = "FALSE" ];
	CALL_SLURM_SCRIPT="rufus_call.slurm"
	SCAFFOLD_SLURM_SCRIPT="rufus_scaffold.slurm"
    echo -e "max_concurrent_jobs=$SLURM_JOB_LIMIT_RUFUS_ARG" >> $CALL_SLURM_SCRIPT
    echo -e "for ((i=0; i<$NUM_CHUNKS; i+=$SLURM_JOB_LIMIT_RUFUS_ARG)); do" >> $CALL_SLURM_SCRIPT
    echo -e "   running_jobs=$(squeue -u $USER | wc -l)" >> $CALL_SLURM_SCRIPT
    echo -e "   while [ "$running_jobs" -ge "$max_concurrent_jobs" ]; do" >> $CALL_SLURM_SCRIPT
    echo -e "       sleep 5m" >> $CALL_SLURM_SCRIPT
    echo -e "       running_jobs=$(squeue -u $USER | wc -l)" >> $CALL_SLURM_SCRIPT
    echo -e "   done" >> $CALL_SLURM_SCRIPT
    echo -e "   step=$(( SLURM_JOB_LIMIT_RUFUS_ARG - 1 ))" >> $CALL_SLURM_SCRIPT
    echo -e "   sbatch --array=${i}-$((i+step)) $SCAFFOLD_SLURM_SCRIPT \$i" >> $CALL_SLURM_SCRIPT
    echo -e "done" >> $CALL_SLURM_SCRIPT
    echo -e "sbatch --depend=afterok:$(squeue -u $USER -o '%A' -h | tr '\n' ':' | sed 's/:$//') $PP_SLURM_SCRIPT" >> $CALL_SLURM_SCRIPT
fi	

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

echo -e "srun singularity exec --bind ${HOST_DATA_DIR_RUFUS_ARG}:/mnt ${CONTAINER_PATH_RUFUS_ARG} bash /opt/RUFUS/post_process/post_process.sh -w $WINDOW_SIZE_RUFUS_ARG -r $REFERENCE_RUFUS_ARG -c $CONTROL_STRING -s $SUBJECT_RUFUS_ARG -d ${BOUND_DATA_DIR}" >> $PP_SLURM_SCRIPT

echo -en "##RUFUS_postProcessCommand=" >> rufus.cmd
echo -e "srun singularity exec --bind ${HOST_DATA_DIR_RUFUS_ARG}:/mnt ${CONTAINER_PATH_RUFUS_ARG} bash /opt/RUFUS/post_process/post_process.sh -w $WINDOW_SIZE_RUFUS_ARG -r $REFERENCE_RUFUS_ARG -c $CONTROL_STRING -s $SUBJECT_RUFUS_ARG -d ${BOUND_DATA_DIR}" >> rufus.cmd
mv rufus.cmd ${HOST_DATA_DIR_RUFUS_ARG}

# Compose invocation script to be executed outside of container
EXE_SCRIPT=launch_rufus.sh
echo -e "#!/bin/bash" > $EXE_SCRIPT
echo -e "" >> $EXE_SCRIPT
echo -e "# This script should be executed after calling the container setup_slurm.sh helper. It requires $PP_SLURM_SCRIPT and $RUFUS_CALL_SCRIPT scripts to be present in the same directory." >> $EXE_SCRIPT 
echo "" >> $EXE_SCRIPT
echo -e "# Launch calling job(s)" >> $EXE_SCRIPT
echo -e "ARRAY_JOB_ID=\$(sbatch --parsable rufus_call.slurm)" >> $EXE_SCRIPT
echo -e "" >> $EXE_SCRIPT
echo -e "# Launch post-process job - will wait on calling phase to complete" >> $EXE_SCRIPT
echo -e "sbatch --depend=afterany:\$ARRAY_JOB_ID $PP_SLURM_SCRIPT" >> $EXE_SCRIPT

# Give user instructions for next step
echo -e "Slurm scripts ready to execute with $EXE_SCRIPT. Please make sure singularity is available in your environment (with e.g. module load singularity), and then run... "
echo -e "bash $EXE_SCRIPT"
