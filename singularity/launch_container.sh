# ONLY SUPPORTS GRCh38 CURRENTLY

# FILL IN VARIABLES
SLURM_ACCT="marth-rw"
SLURM_PARTITION="marth-shared-rw"
WORKING_DIR=$(pwd)
EMAIL="u0746015@utah.edu"
JOB_NAME="RUFUS_TEST"
TIME_LIMIT="3-00:00:00"
SUBJECT_FILE=data/NA12878-NA18517-5percent-lib1-umi-liquid_tumor.sorted.bam
CONTROL_FILE=data/NA12878.sorted.bam
KMER_CUTOFF=5
CHUNK_SIZE=1000000 # set to 1 for single, complete run
NUM_NODES_PER_JOB=1
JOB_LIMIT=900

# LEFTOFF HERE
# todo: put this directly in here - hard coded to GRCh38
NUM_CHUNKS=$(bash utilities/get_num_chunks.sh $CHUNK_SIZE)


cd $WORKING_DIR
mkdir slurm_out
mkdir -p slurm_out/individual_jobs

header_lines=("#!/bin/bash"
"#SBATCH --time=${TIME_LIMIT}" 
"#SBATCH --account=${SLURM_ACCT}" 
"#SBATCH --partition=${SLURM_PARTITION}"
"#SBATCH --nodes=${NUM_NODES_PER_JOB}"
"#SBATCH -o ${WORKING_DIR}/slurm_out/%j.out"
"#SBATCH -e ${WORKING_DIR}/slurm_err/%j.err"
"#SBATCH --mail-type=ALL" 
"#SBATCH --mail-user=u0746015@utah.edu"
"#SBATCH -a 0-${NUM_CHUNKS}%${JOB_LIMIT}"
)

SLURM_SCRIPT="run_rufus.slurm"
for line in "${header_lines[@]}"
do
    echo -e "$line" >> $SLURM_SCRIPT
done
echo "\n" >> $SLURM_SCRIPT

# todo: need to test getChunk.sh
echo "REGION_ARG=$(bash getChunk.sh ${SLURM_ARRAY_TASK_ID})" >> $SLURM_SCRIPT

# todo: need to edit this for singularity setup
echo -e "bash \$RUFUS_ROOT/runRufus.sh -s $SUBJECT -c $CONTROL -r $REFERENCE \
 -f ${RESOURCE_DIR}GRCh38_full_analysis_set_plus_decoy_hla.25.Jhash -m $KMER_CUTOFF \
  -k 25 -t 40 -L -vs $REGION_ARG -v $SLURM_ARRAY_TASK_ID" >> $SLURM_SCRIPT

sbatch $SLURM_SCRIPT

n=${#CHRS[@]}
for (( i = 0; i < n; i++ ))
do
    curr_len="${CHR_LENGTHS[i]}"
    curr_chr="${CHRS[i]}"
    ((chunks=curr_len/CHUNK_SIZE+1))
    echo "number of chunks is $chunks"
 
    start_coord=1
    echo "Starting loop for index $i chr${curr_chr}"
    for (( c = 1; c <= chunks; c++ ))
    do
	    curr_jobs=($(squeue -u $UNID | wc -l))
        while [ $curr_jobs -gt 990 ]; do
		    sleep 60;
	        curr_jobs=($(squeue -u $UNID | wc -l))
	    done

        ((end_coord=start_coord+CHUNK_SIZE-1))
        if [ $end_coord -gt $curr_len ]; then
            echo "passed end coordinate, setting to $curr_len"
            end_coord=$curr_len
        fi
	
	    # Make region dir
        TARGET_DIR="chr${curr_chr}/chr${curr_chr}_${start_coord}_${end_coord}/"
        mkdir -p $TARGET_DIR
        
	    # Compose slurm script
	    SLURM_FILE="${TARGET_DIR}rufus_HCC.sh"
	    echo -e "$SLURM_SCRIPT" > $SLURM_FILE 
        echo -e "cd ${OUT_DIR}${TARGET_DIR} \n" >> $SLURM_FILE
	    echo -e "bash \$RUFUS_ROOT/runRufus.sh -s $SUBJECT -c $CONTROL -r $REFERENCE -f ${RESOURCE_DIR}GRCh38_full_analysis_set_plus_decoy_hla.25.Jhash -m $KMER_CUTOFF -k 25 -t 40 -L -vs -R chr${curr_chr}:${start_coord}-${end_coord}" >> $SLURM_FILE 
        
        # Get run going
	    sbatch $SLURM_FILE

        # Advance start coordinate
        ((start_coord=end_coord+1))
    done
    echo "Done with chr${curr_chr}"
done

# todo: how do I ensure that wait for all jobs to complete here?
#todo: what args do I need to pass here?
sbatch --depend=afterany:$SLURM_ARRAY_JOB_ID singularity exec bash postProcessChunkRun.sh 
