# ONLY SUPPORTS GRCh38 CURRENTLY

#todo: put all user fill in into a .args file
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
CHUNK_SIZE=1000000 # set to 0 for whole genome

# Making into variables for now so can update later as a function of arg genome build
GENOME_BUILD="GRCh38"
ref_hash="GRCh38_full_analysis_set_plus_decoy_hla.25.Jhash"

# FILL IN SAMPLE VARIABLES
SUBJECT_FILE=data/NA12878-NA18517-5percent-lib1-umi-liquid_tumor.sorted.bam
CONTROL_FILE=data/NA12878.sorted.bam
RESOURCE_DIR="/resources"

# Load resource functions and arguments
#source rufus.args - TODO: this will need to change to singularity path
. /Users/genetics/Documents/code/RUFUS/singularity/launch_utilities/chunk_utilities.sh

num_chunks=$(get_num_chunks "$CHUNK_SIZE" "$GENOME_BUILD")
echo "num chunks is $num_chunks"

WORKING_DIR=$(pwd)
cd "$WORKING_DIR" || exit
mkdir -p slurm_out

# Compose header
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
)
if ((num_chunks > 1)); then
    header_lines+=("#SBATCH -a 0-${num_chunks}%${JOB_LIMIT}")
fi

# Compose individual SLURM scripts
SLURM_SCRIPT="run_rufus.slurm"
for line in "${header_lines[@]}"
do
    echo -e "$line" >> $SLURM_SCRIPT
done
echo "\n" >> $SLURM_SCRIPT

# Add in region argument if there are multiple chunks
if ((num_chunks > 1)); then
    echo "source ./chunk_utilites.sh" >> $SLURM_SCRIPT
    echo "region_arg=$(get_chunk_region "$SLURM_ARRAY_TASK_ID" "$CHUNK_SIZE" "$GENOME_BUILD")" >> $SLURM_SCRIPT
    echo "REGION_ARG= -R ${region_arg}" >> $SLURM_SCRIPT
else
    echo "REGION_ARG=\"\" >> $SLURM_SCRIPT"
fi

# Add in the rufus command
echo -en "bash \$RUFUS_ROOT/runRufus.sh -s $SUBJECT_FILE -c $CONTROL_FILE -r $REFERENCE \
 -f ${RESOURCE_DIR}${ref_hash} -m $KMER_CUTOFF \
  -k 25 -t 40 -L -vs $REGION_ARG" >> $SLURM_SCRIPT
echo -e "-v $SLURM_ARRAY_TASK_ID" >> $SLURM_SCRIPT

# Submit the SLURM script
# todo: how do I submit a slurm array
#sbatch $SLURM_SCRIPT

# Wait here until all scripts done
#todo: what args do I need to pass here?
#sbatch --depend=afterany:$SLURM_ARRAY_JOB_ID singularity exec bash postProcessChunkRun.sh
