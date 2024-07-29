#!/bin/bash

apt install -y nfs-common
echo "installed all tools and dependencies"

# mount nfs drive
mount 172.31.17.24:/mnt /mnt
echo "mounted nfs drive"

# recompile RUFUS to fix jellyfish library link
cd ~/RUFUS
rm -rf bin
mkdir bin
cd bin
cmake ..
make

# check to make sure paths are updated
# DO THIS DON'T FORGET AGAIN
echo "CHECK generate_mb_scripts_ec2.sh to make sure looking for resources in /mnt/resources now..."
echo "UPDATE RUFUS KMER PARAMETER then run the rest of setup script"

# run generator
cd 
./generate_mb_scripts_ec2.sh /mnt/SMHTCOLO829BLT50-X-X-M45-A001-dac-SMAARNRVZGBE-insilico500X_GRCh38.aligned.sorted.bam /mnt/COLO829BL_Ill_230X.bam 

# run RUFUS
cd /mnt/rufus-runs/SMHTCOLO829BLT50-X-X-M45-A001-dac-SMAARNRVZGBE-insilico500X_GRCh38.aligned.sorted/
mkdir logs
cd launchers
ls | parallel -j 40 'bash {} >/mnt/rufus-runs/SMHTCOLO829BLT50-X-X-M45-A001-dac-SMAARNRVZGBE-insilico500X_GRCh38.aligned.sorted_k3/logs/{}.out 2>/mnt/rufus-runs/SMHTCOLO829BLT50-X-X-M45-A001-dac-SMAARNRVZGBE-insilico500X_GRCh38.aligned.sorted_k3/logs/{}.err'
