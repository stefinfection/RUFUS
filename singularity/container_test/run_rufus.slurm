#!/bin/bash
#SBATCH --job-name=RUFUS
#SBATCH --time=01:00:00
#SBATCH --account=marth-rw
#SBATCH --partition=marth-shared-rw
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH -o /home/ubuntu/container_test/slurm_out/%A_%a.out
#SBATCH -e /home/ubuntu/container_test/slurm_err/%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=u0746015@utah.edu
#SBATCH -a 0-3102%20
#SBATCH --nodes=1-20
#SBATCH --mem=64G

region_arg=$(singularity exec rufus.sif bash /opt/RUFUS/singularity/launch_utilities/get_region.sh "$SLURM_ARRAY_TASK_ID" "1000" "GRCh38")
REGION_ARG="-R $region_arg"
srun singularity exec --bind /home/ubuntu/RUFUS/testRun//:/mnt rufus.sif bash /opt/RUFUS/runRufus.sh -s /mnt/Child.bam -c Mother.bam -c Father.bam -f /opt/RUFUS/resources/references/prebuilt_hashes/GRCh38_full_analysis_set_plus_decoy_hla.25.Jhash -r fake_ref.fa -m 12 -k 25 -t 40 -L -vs $REGION_ARG
