#!/bin/bash
#SBATCH --job-name=RUFUS_TEST
#SBATCH --time=3-00:00:00
#SBATCH --account=marth-rw
#SBATCH --partition=marth-shared-rw
#SBATCH --nodes=1
#SBATCH -o /Users/genetics/Documents/code/RUFUS/singularity/launch_utilities/slurm_out/%j.out
#SBATCH -e /Users/genetics/Documents/code/RUFUS/singularity/launch_utilities/slurm_err/%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=u0746015@utah.edu
#SBATCH -a 0-3102%900

. /Users/genetics/Documents/code/RUFUS/singularity/launch_utilities/chunk_utilities.sh

region_arg=$(get_chunk_region $SLURM_ARRAY_TASK_ID "1000000" "GRCh38")
REGION_ARG="-R $region_arg"
bash $RUFUS_ROOT/runRufus.sh -s data/NA12878-NA18517-5percent-lib1-umi-liquid_tumor.sorted.bam -c data/NA12878.sorted.bam -r /scratch/ucgd/lustre/work/marth/shared/references/human/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa -f /resources/GRCh38_full_analysis_set_plus_decoy_hla.25.Jhash -m 5 -k 25 -t 40 -L -vs $REGION_ARG -v $SLURM_ARRAY_TASK_ID
