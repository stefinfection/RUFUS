if [ -e /opt/RUFUS/runRufus.sh ]
then 
	echo "starting test"
else
	echo "could not find /opt/runRufus.sh, are you in the RUFUS/runTest/ directory?"
	exit
fi 

RUFUS_PATH=/opt/RUFUS
RUFUS_TEST_PATH=/mnt
${RUFUS_PATH}/runRufus.sh -s ${RUFUS_TEST_PATH}/Child.bam -c ${RUFUS_TEST_PATH}/Mother.bam -c ${RUFUS_TEST_PATH}/Father.bam -k 25 -t 40  -r ${RUFUS_PATH}/resources/references/small_test_human_reference_v37_decoys.fa $1
