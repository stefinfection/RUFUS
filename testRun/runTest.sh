if [ -e /home/ubuntu/RUFUS/runRufus.sh ]
then 
	echo "starting test"
else
	echo "could not find ./../runRufus.sh, are you in the RUFUS/runTest/ directory?"
	exit
fi 

RUFUS_PATH=/home/ubuntu/RUFUS
RUFUS_TEST_PATH="${RUFUS_PATH}/testRun"
${RUFUS_PATH}/runRufus.sh -s ${RUFUS_TEST_PATH}/Child.bam -c ${RUFUS_TEST_PATH}/Mother.bam -c ${RUFUS_TEST_PATH}/Father.bam -k 25 -t 40  -r ${RUFUS_PATH}/resources/references/small_test_human_reference_v37_decoys.fa $1
