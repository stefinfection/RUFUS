#!/bin/bash
bash /home/ubuntu/RUFUS/singularity/launch_container.sh -d /home/ubuntu/RUFUS/testRun/ -s Child.bam -c Mother.bam,Father.bam -r fake_ref.fa -a marth-rw -p marth-shared-rw -w 1000 -t "00:05:00" -m 12 -l 20 -z 40 -e "u0746015@utah.edu"
