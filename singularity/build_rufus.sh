export SINGULARITY_TMPDIR=/home/ubuntu/singularity_tmp
export SINGULARITY_BIND=/home/ubuntu/rufus_data:/mnt
<<<<<<< HEAD
sudo -E singularity build --sandbox rufus.sif rufus.def
=======
sudo -E singularity build rufus.sif rufus.def
>>>>>>> 927c60278967e8c3940f191d74f3b7fe96aefa0c
