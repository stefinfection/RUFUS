DATA_PATH_HOST=/home/ubuntu/RUFUS/testRun

#todo: this will be programmatically generated based on rufus.config

singularity exec --bind $DATA_PATH_HOST:/mnt /home/ubuntu/rufus.sif bash /opt/RUFUS/runRufus.sh -s /mnt/Child.bam -c /mnt/Mother.bam -c /mnt/Father.bam -k 25 -t 40  -r /opt/RUFUS/resources/references/small_test_human_reference_v37_decoys.fa
