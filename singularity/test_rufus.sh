RUFUS_PATH=/home/ubuntu/RUFUS
RUFUS_TEST_PATH="${RUFUS_PATH}/testRun"
singularity exec /home/ubuntu/rufus.sif bash ${RUFUS_PATH}/runRufus.sh -s ${RUFUS_TEST_PATH}/Child.bam -c ${RUFUS_TEST_PATH}/Mother.bam -c ${RUFUS_TEST_PATH}/Father.bam -k 25 -t 40  -r ${RUFUS_PATH}/resources/references/small_test_human_reference_v37_decoys.fa
