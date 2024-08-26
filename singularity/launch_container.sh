# ONLY SUPPORTS GRCh38 CURRENTLY

# FILL IN SLURM VARIABLES
SLURM_ACCT="marth-rw"
SLURM_PARTITION="marth-shared-rw"
EMAIL="u0746015@utah.edu"
JOB_NAME="RUFUS_TEST"
TIME_LIMIT="3-00:00:00"
NUM_NODES_PER_JOB=1
JOB_LIMIT=900

# FILL IN RUFUS VARIABLES
KMER_CUTOFF=5
CHUNK_SIZE=1000000 # set to 1 for single, complete run

# FILL IN SAMPLE VARIABLES
SUBJECT_FILE=data/NA12878-NA18517-5percent-lib1-umi-liquid_tumor.sorted.bam
CONTROL_FILE=data/NA12878.sorted.bam

# LEFTOFF HERE
# todo: put this directly in here - hard coded to GRCh38
NUM_CHUNKS=$(bash utilities/get_num_chunks.sh $CHUNK_SIZE)

WORKING_DIR=$(pwd)
cd "$WORKING_DIR" || exit
mkdir slurm_out
mkdir -p slurm_out/individual_jobs

header_lines=("#!/bin/bash"
"#SBATCH --job-name=${JOB_NAME}"
"#SBATCH --time=${TIME_LIMIT}" 
"#SBATCH --account=${SLURM_ACCT}" 
"#SBATCH --partition=${SLURM_PARTITION}"
"#SBATCH --nodes=${NUM_NODES_PER_JOB}"
"#SBATCH -o ${WORKING_DIR}/slurm_out/%j.out"
"#SBATCH -e ${WORKING_DIR}/slurm_err/%j.err"
"#SBATCH --mail-type=ALL" 
"#SBATCH --mail-user=${EMAIL}"
"#SBATCH -a 0-${NUM_CHUNKS}%${JOB_LIMIT}"
)

SLURM_SCRIPT="run_rufus.slurm"
for line in "${header_lines[@]}"
do
    echo -e "$line" >> $SLURM_SCRIPT
done
echo "\n" >> $SLURM_SCRIPT

# todo: need to test getChunk.sh
echo "REGION_ARG=$(bash getChunk.sh "${SLURM_ARRAY_TASK_ID}")" >> $SLURM_SCRIPT

# todo: need to edit this for singularity setup
echo -e "bash \$RUFUS_ROOT/runRufus.sh -s $SUBJECT_FILE -c $CONTROL_FILE -r $REFERENCE \
 -f ${RESOURCE_DIR}GRCh38_full_analysis_set_plus_decoy_hla.25.Jhash -m $KMER_CUTOFF \
  -k 25 -t 40 -L -vs $REGION_ARG -v $SLURM_ARRAY_TASK_ID" >> $SLURM_SCRIPT

sbatch $SLURM_SCRIPT

#todo: what args do I need to pass here?
sbatch --depend=afterany:$SLURM_ARRAY_JOB_ID singularity exec bash postProcessChunkRun.sh 
