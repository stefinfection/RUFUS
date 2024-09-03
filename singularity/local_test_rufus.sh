#!/bin/bash
DATA_PATH_LOCAL=/home/ubuntu/RUFUS/testRun
RUFUS_PATH_LOCAL=/home/ubuntu/RUFUS

bash /home/ubuntu/RUFUS/runRufus.sh -s ${DATA_PATH_LOCAL}/Child.bam -c ${DATA_PATH_LOCAL}/Mother.bam -c ${DATA_PATH_LOCAL}/Father.bam -k 25 -t 40  -r ${RUFUS_PATH_LOCAL}/resources/references/small_test_human_reference_v37_decoys.fa -local
