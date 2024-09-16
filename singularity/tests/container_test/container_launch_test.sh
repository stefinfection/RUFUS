#!/bin/bash
PATH_TO_RUFUS_CONTAINER=/home/ubuntu/
DATA_PATH_HOST=/home/ubuntu/RUFUS/testRun/

singularity exec --bind $DATA_PATH_HOST:/mnt $PATH_TO_RUFUS_CONTAINER/rufus.sif bash /opt/RUFUS/singularity/launch_container.sh -d /mnt/ -s Child.bam -c Mother.bam -c Father.bam -r fake_ref.fa 
